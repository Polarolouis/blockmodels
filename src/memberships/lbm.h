
struct LBM
{
    mat Z1;
    mat Z2;
    rowvec alpha1;
    rowvec alpha2;

    class network
    {
    public:
        /* Here you should put all the variable which define the network
         * for scalar model, a mat field is good choice.
         * See bernoulli, poisson model for example.
         *
         * Covariates must also stocked here.
         * See poisson_covariates model for example
         */

        mat adj;

        bool has_edge_covariates{false};
        cube covariates;

        bool has_row_covariates{false};
        bool has_col_covariates{false};
        mat row_nodes_covariates;
        mat col_nodes_covariates;

        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */
        mat adjZD;
        mat adjZDt;
        mat MonesZD;

        mat Mones;
        mat Monest;
        mat adjt;

        network(Rcpp::List &network_from_R)
        {
            /* Here you must define how initialize the fields describing the
             * network. Precomputed value must be computed here.
             * See bernoulli, poisson or poisson_covariates for example.
             *
             * For scalar network, the provided list have a adjacency field
             * which contains the adjacency matrix.
             * For scalar network with covariates vectors on edges (existent or
             * not), the provided list have a covariates field which contain a
             * list of matrices. This matrices have the same size than the
             * adjacency matrix, and the i-th matrix is the matrix of the i-th
             * covariate on all edges.
             */
            adj = Rcpp::as<mat>(network_from_R["adjacency"]);

            if (network_from_R["nodes_covariates"].containsElementNamed("row"))
            {
                are_row_covariates = true;
                row_nodes_covariates = Rcpp::as<mat>(network_from_R["nodes_covariates"]["row"]);
                row_nodes_covariates.set_size()
            }
            if (network_from_R["nodes_covariates"].containsElementNamed("col"))
            {
                are_col_covariates = true;
                col_nodes_covariates = Rcpp::as<mat>(network_from_R["nodes_covariates"]["col"]);
            }

            if (network_from_R.containsElementNamed("covariates"))
            {
                Rcpp::List covariates_list = network_from_R["covariates"];

                covariates.set_size(adj.n_rows, adj.n_cols, covariates_list.size());
                for (int k = 0; k < covariates_list.size(); k++)
                    covariates.slice(k) = Rcpp::as<mat>(covariates_list[k]);
            }
            adjZD = fill_diag(adj, 0);
            Mones = ones<mat>(adj.n_rows, adj.n_cols);
            MonesZD = fill_diag(Mones, 0);
            Monest = Mones.t();
        }
    };

    LBM(Rcpp::List &membership_from_R)
    {
        mat origZ1 = membership_from_R["Z1"];
        mat origZ2 = membership_from_R["Z2"];
        Z1 = origZ1;
        Z2 = origZ2;
        double tol1 = TOL_P / Z1.n_rows;
        double tol2 = TOL_P / Z2.n_rows;
        boundaries(Z1, tol1, 1 - tol1);
        boundaries(Z2, tol2, 1 - tol2);
        Z1 /= repmat(sum(Z1, 1), 1, Z1.n_cols);
        Z2 /= repmat(sum(Z2, 1), 1, Z2.n_cols);
        alpha1 = sum(Z1, 0) / Z1.n_rows;
        alpha2 = sum(Z2, 0) / Z2.n_rows;
    }

    LBM &operator=(const LBM &orig)
    {
        Z1 = orig.Z1;
        Z2 = orig.Z2;
        alpha1 = orig.alpha1;
        alpha2 = orig.alpha2;

        return *this;
    }

    inline double entropy()
    {
        return -accu(Z1 % log(Z1)) - accu(Z2 % log(Z2));
    }

    template <class model_type, class network_type>
    inline void e_step(model_type &model, network_type &net)
    {
        double tol1 = TOL_P / Z1.n_rows;
        double tol2 = TOL_P / Z2.n_rows;
        double step_size;
        unsigned int niter = 0;

        do
        {
            mat lZ1 = repmat(log(alpha1), Z1.n_rows, 1);
            mat lZ2 = repmat(log(alpha2), Z2.n_rows, 1);

            e_fixed_step(*this, model, net, lZ1, lZ2);

            lZ1 -= repmat(mean(lZ1, 1), 1, lZ1.n_cols);
            lZ2 -= repmat(mean(lZ2, 1), 1, lZ2.n_cols);

            lZ1 -= repmat(max(lZ1, 1), 1, lZ1.n_cols);
            lZ2 -= repmat(max(lZ2, 1), 1, lZ2.n_cols);

            lZ1 = exp(lZ1);
            lZ2 = exp(lZ2);

            lZ1 /= repmat(sum(lZ1, 1), 1, lZ1.n_cols);
            lZ2 /= repmat(sum(lZ2, 1), 1, lZ2.n_cols);

            boundaries(lZ1, tol1, 1 - tol1);
            boundaries(lZ2, tol2, 1 - tol2);
            lZ1 /= repmat(sum(lZ1, 1), 1, lZ1.n_cols);
            lZ2 /= repmat(sum(lZ2, 1), 1, lZ2.n_cols);

            double step_size1 = max(max(abs(Z1 - lZ1)));
            double step_size2 = max(max(abs(Z2 - lZ2)));

            step_size = (step_size1 > step_size2) ? step_size1 : step_size2;

#ifdef DEBUG_E
            fprintf(stderr, "E iteration: %i %f [%f %f]\n", niter, step_size, step_size1, step_size2);
#endif

            niter++;

            Z1 = lZ1;
            Z2 = lZ2;
        } while (step_size > TOL_F && niter < 10);
    }

    inline double m_step()
    {
        alpha1 = sum(Z1, 0) / Z1.n_rows;
        alpha2 = sum(Z2, 0) / Z2.n_rows;

        return accu(Z1 * log(alpha1).t()) + accu(Z2 * log(alpha2).t());
    }

    inline Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["Z1"] = Z1;
        values["alpha1"] = alpha1;
        values["Z2"] = Z2;
        values["alpha2"] = alpha2;

        return values;
    }
};
