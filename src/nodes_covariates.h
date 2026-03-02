inline mat softmax(const mat &Z)
{
    mat Zshift = Z;
    vec maxZ = max(Z, 1); // row max

    Zshift.each_col() -= maxZ; // stability

    mat P = exp(Zshift);
    vec rowsum = sum(P, 1);

    P.each_col() /= rowsum;

    return P;
}

inline mat log_softmax(const mat &Z)
{
    vec maxZ = max(Z, 1);
    mat Zshift = Z;
    Zshift.each_col() -= maxZ;

    vec logsumexp = log(sum(exp(Zshift), 1)) + maxZ;

    mat logP = Z;
    logP.each_col() -= logsumexp;

    return logP;
}

inline mat build_full_B(const mat &Btilde)
{
    mat B(Btilde.n_rows, Btilde.n_cols + 1, fill::zeros);
    B.cols(0, Btilde.n_cols - 1) = Btilde;
    return B;
}

double objective_gradient(
    const mat &X,
    const mat &T,
    mat &Btilde,
    mat &grad)
{
    mat B = build_full_B(Btilde);

    mat Z = X * B;

    mat logP = log_softmax(Z);
    double L = accu(T % logP);

    mat P = softmax(Z);

    vec s = sum(T, 1);
    mat SP = P;
    SP.each_col() %= s;

    mat Gfull = X.t() * (T - SP);

    // keep only identifiable parameters
    grad = Gfull.cols(0, Btilde.n_cols - 1);

    return L;
}
double line_search_bfgs(
    const mat &X,
    const mat &T,
    const mat &Btilde,
    const mat &grad,
    const mat &direction,
    double L)
{
    double alpha = 1.0;
    const double c1 = 1e-4;
    const double slope = accu(grad % direction);

    if (slope <= 0.0)
        return 0.0;

    for (int i = 0; i < 25; i++)
    {
        mat Bnew = Btilde + alpha * direction;
        mat gtmp;
        double Lnew = objective_gradient(X, T, Bnew, gtmp);

        // Armijo condition for maximization
        if (Lnew >= L + c1 * alpha * slope)
            return alpha;

        alpha *= 0.5;
    }

    return 0.0;
}

mat optimize_softmax(
    const mat &X,
    const mat &T,
    int max_iter = 1000,
    double tol = 1e-4)
{
    int p = X.n_cols;
    int R = T.n_cols;
    int q = R - 1;
    int nparam = p * q;
    int last_iter = 0;
    // only R-1 columns optimized
    mat Btilde(p, q, fill::zeros);
    mat grad;

    double L = objective_gradient(X, T, Btilde, grad);

    mat H = eye<mat>(nparam, nparam); // inverse Hessian approximation

    for (int iter = 0; iter < max_iter; iter++)
    {
        vec g = vectorise(grad);

        if (norm(g, 2) < tol)
            break;

        vec d = H * g; // ascent direction (maximize)

        if (dot(g, d) <= 0.0)
        {
            d = g;
            H.eye();
        }

        mat direction = reshape(d, p, q);
        double step = line_search_bfgs(X, T, Btilde, grad, direction, L);

        if (step <= 0.0)
            break;

        mat Bold = Btilde;
        vec gold = g;

        Btilde += step * direction;

        double Lnew = objective_gradient(X, T, Btilde, grad);
        vec gnew = vectorise(grad);

        double rel_obj_diff = std::abs(Lnew - L) / (1.0 + std::abs(L));

        vec s = vectorise(Btilde - Bold);
        vec y = gnew - gold;
        double ys = dot(y, s);

        if (ys > 1e-12)
        {
            double rho = 1.0 / ys;
            mat I = eye<mat>(nparam, nparam);
            mat syT = s * y.t();
            mat ysT = y * s.t();
            H = (I - rho * syT) * H * (I - rho * ysT) + rho * (s * s.t());
        }
        else
        {
            H.eye();
        }

        L = Lnew;
        last_iter = iter;

        if (rel_obj_diff < tol)
            break;

        if (norm(gnew, 2) < tol)
            break;
        #ifdef DEBUG_M
        if (iter % 50 == 0 || iter == max_iter - 1) {
        Rcpp::Rcout << "Iteration " << iter + 1 << "/" << max_iter << ": L = " << Lnew << endl;
        }
        #endif
    }

    if (last_iter >= max_iter) {
        Rcpp::warning("Optimization for nodes covariates did not converge after %i steps.", last_iter);
    }

    return build_full_B(Btilde);
}