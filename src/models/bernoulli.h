
class bernoulli
{
    public:

    class network
    {
        public:
        mat adj; // adjacency matrix
        mat maskNA;

        // precalculated matrices for SBM
        mat adjZD; 
        mat ones_minus_adj_ZD;
        mat adjZDt;
        mat ones_minus_adj_ZDt;
        mat onesZD;

        // precalculated matrices for LBM
        mat ones_minus_adj;
        mat adjt;
        mat ones_minus_adjt;
        mat adj_ones;


        network(Rcpp::List & network_from_R)
        {
            mat adj_orig = network_from_R["adjacency"];
            double na_replace_value = 0;
            if(network_from_R.containsElementNamed("na_replace_value"))
            {
                na_replace_value = Rcpp::as<double>(network_from_R["na_replace_value"]);
            }

            maskNA = compute_mask(adj_orig);
            adj = replace_missing_values(adj_orig, na_replace_value) % maskNA;
            adjZD = fill_diag(adj,0);
            ones_minus_adj_ZD = fill_diag(maskNA-adj,0);
            adjZDt = adjZD.t();
            ones_minus_adj_ZDt = ones_minus_adj_ZD.t();
            onesZD = fill_diag(maskNA,0);

            ones_minus_adj = maskNA-adj;
            adjt = adj.t();
            ones_minus_adjt = ones_minus_adj.t();
            adj_ones = maskNA;

        }
    };

    /* parameters */
    unsigned int n_parameters;
    mat pi;

    bernoulli(SBM & membership, bernoulli::network & net)
    {
        n_parameters = membership.Z.n_cols * membership.Z.n_cols;
        pi.set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    bernoulli(SBM_sym & membership, bernoulli::network & net)
    {
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols+1) / 2;
        pi.set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    bernoulli(LBM & membership, bernoulli::network & net)
    {
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols;
        pi.set_size(membership.Z1.n_cols,membership.Z2.n_cols);
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["pi"] = pi;
        values["n_parameters"] = n_parameters;

        return values;
    }
};

template<>
inline
void e_fixed_step(SBM & membership, bernoulli & model, bernoulli::network & net, mat & lZ)
{
    lZ+= net.adjZD * membership.Z * log(model.pi).t()
       + net.ones_minus_adj_ZD * membership.Z * log(1-model.pi).t()
       + net.adjZDt * membership.Z * log(model.pi)
       + net.ones_minus_adj_ZDt * membership.Z * log(1-model.pi);
}

template<>
inline
void e_fixed_step(SBM_sym & membership, bernoulli & model, bernoulli::network & net, mat & lZ)
{
    lZ+= net.adjZD * membership.Z * log(model.pi)
       + net.ones_minus_adj_ZD * membership.Z * log(1-model.pi);
}

template<>
inline
void e_fixed_step(LBM & membership, bernoulli & model, bernoulli::network & net, mat & lZ1, mat & lZ2)
{
    lZ1 += net.adj * membership.Z2 * log(model.pi).t()
         + net.ones_minus_adj * membership.Z2 * log(1-model.pi).t();
    lZ2 += net.adjt * membership.Z1 * log(model.pi)
         + net.ones_minus_adjt * membership.Z1 * log(1-model.pi);
}

template<>
inline
double m_step(SBM & membership, bernoulli & model, bernoulli::network & net)
{
    mat effective_counts = membership.Z.t() * net.onesZD * membership.Z;
    if(accu(effective_counts)<=0)
    {
        Rcpp::stop("No valid non-missing off-diagonal adjacency values in bernoulli SBM network.");
    }
    mat safe_counts = effective_counts;
    safe_counts.elem(find(safe_counts<=0)).ones();

    model.pi = (membership.Z.t() * net.adjZD * membership.Z)
                /
               safe_counts;

    return
        (
            accu(
                    (log(model.pi)-log(1-model.pi)) 
                    % 
                    (membership.Z.t() * net.adjZD * membership.Z)
                )
            +
            accu(
                    log(1-model.pi)
                    %
                        effective_counts
                )
        );
}

template<>
inline
double m_step(SBM_sym & membership, bernoulli & model, bernoulli::network & net)
{
    return m_step<SBM>(membership,model,net)/2;
}

template<>
inline
double m_step(LBM & membership, bernoulli & model, bernoulli::network & net)
{
    mat effective_counts = membership.Z1.t() * net.adj_ones * membership.Z2;
    if(accu(effective_counts)<=0)
    {
        Rcpp::stop("No valid non-missing adjacency values in bernoulli LBM network.");
    }
    mat safe_counts = effective_counts;
    safe_counts.elem(find(safe_counts<=0)).ones();

    model.pi = (membership.Z1.t() * net.adj * membership.Z2)
                /
               safe_counts;

    return
        (
            accu(
                    log(model.pi) 
                    % 
                    (membership.Z1.t() * net.adj * membership.Z2)
                )
            +
            accu(
                    log(1-model.pi)
                    %
                        (effective_counts - (membership.Z1.t() * net.adj * membership.Z2))
                )
        );
}


