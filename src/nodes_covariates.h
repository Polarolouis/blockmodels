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
double line_search(
    const mat &X,
    const mat &T,
    mat &Btilde,
    const mat &grad,
    double L)
{
    double alpha = 1.0, c = 1e-4;

    for (int i = 0; i < 20; i++)
    {
        mat Bnew = Btilde + alpha * grad;
        mat gtmp;
        double Lnew = objective_gradient(X, T, Bnew, gtmp);

        if (Lnew >= L + c * alpha * accu(grad % grad))
            return alpha;

        alpha *= 0.5;
    }
    return alpha;
}

mat optimize_softmax(
    const mat &X,
    const mat &T,
    int max_iter = 500,
    double tol = 1e-6)
{
    int p = X.n_cols;
    int R = T.n_cols;

    // only R-1 columns optimized
    mat Btilde(p, R - 1, fill::zeros);
    mat grad;

    double L = objective_gradient(X, T, Btilde, grad);

    for (int iter = 0; iter < max_iter; iter++)
    {
        double step = line_search(X, T, Btilde, grad, L);

        Btilde += step * grad;

        double Lnew = objective_gradient(X, T, Btilde, grad);

        if (norm(grad, "fro") < tol)
            break;

        L = Lnew;
    }

    return build_full_B(Btilde);
}