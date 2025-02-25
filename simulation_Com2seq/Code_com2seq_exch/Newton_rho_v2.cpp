#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube CalculateXY(const arma::mat & X, const arma::mat & Y){
    int K = X.n_cols;
    int n_sam = X.n_rows;
    arma::cube XY(K, K, n_sam);
    
    for (int s = 0; s < n_sam; s++){
        for (int i = 0; i < K; i++){
            XY(i, i, s) = X(s, i) * Y(s, i);
            for (int j = i+1; j < K; j++){
                XY(i, j, s) = X(s, i) * Y(s, j);
                XY(j, i, s) = XY(i, j, s);
            }
        }
    }
    return XY;
}//outer product of X and Y(e.g., X and (-1,1,\tilde{0}), X and X)

// estimate rho
arma::vec compute_residual(const arma::vec& freq_table, const arma::vec& u, const arma::vec& w) {
    
    return sqrt(w) % (freq_table - u) / sqrt(u % (1 - u));
}

// [[Rcpp::export]]
List Newton_cor(arma::mat freq_table, arma::mat N_denominator, arma::mat X, 
                arma::mat beta_init, arma::mat weight,
                double tol, int iter_max, double Firth_thresh, bool robust_var,
                arma::vec prop_presence, bool get_var) {

    int n_otu = freq_table.n_cols;
    int n_sam = freq_table.n_rows;
    int n_sam_1 = n_sam/2;
    int K = X.n_cols;
    
    arma::mat X1 = X.submat(0,0,n_sam_1-1,K-1);
    arma::mat X2 = X.submat(n_sam_1,0,n_sam-1,K-1);
    arma::cube X1X1  = CalculateXY(X1, X1);
    arma::cube X2X2  = CalculateXY(X2, X2);
    arma::cube X1X2  = CalculateXY(X1, X2);
    arma::cube X2X1  = CalculateXY(X2, X1);
    arma::mat beta = beta_init;
    arma::vec z;
    arma::vec u;
    arma::vec v;
    arma::vec v1;
    arma::vec v2;
    arma::vec u1;
    arma::vec u2;
    arma::vec DV; 
    arma::vec DV1; 
    arma::vec DV2; 
    arma::vec S;
    arma::vec J_temp11;
    arma::vec J_temp21;
    arma::vec J_temp22;
    arma::vec J_temp12;
    arma::vec J_temp11_1;
    arma::vec J_temp21_1;
    arma::vec J_temp22_1;
    arma::vec J_temp12_1;
    arma::vec H_temp11;
    arma::vec H_temp21;
    arma::vec H_temp22;
    arma::vec H_temp12;
    arma::mat J(K, K, arma::fill::zeros);
    arma::mat J_inv;
    arma::mat Sigma(K, K, arma::fill::zeros);
    arma::vec H_temp;
    arma::mat H(K, K, arma::fill::zeros);
    arma::vec step;
    arma::cube variance;
    arma::vec rho_estimates(n_otu);
    arma::vec final_step_sum(n_otu);

    if (get_var) variance = arma::cube(K, K, n_otu, arma::fill::zeros);
    
    for (int j = 0; j < n_otu; j++){

        beta.col(j) = beta_init.col(j);

        for (int i = 0; i < iter_max; i++){

            z = exp(X * beta.col(j));
            u = z / (1 + z);    
            
            v =  N_denominator.col(j) % u % (1 - u);// for the jth taxon
            arma::vec v1 = v.subvec(0, n_sam_1 - 1);
            arma::vec v2 = v.subvec(n_sam_1, n_sam - 1); 
            arma::vec u1 = u.subvec(0, n_sam_1 - 1);
            arma::vec u2 = u.subvec(n_sam_1, n_sam - 1); 
            // estimate rho
            
            arma::vec freq_table1_j = freq_table.submat(0, j, n_sam_1 - 1, j);
            arma::vec freq_table2_j = freq_table.submat(n_sam_1, j, n_sam - 1, j);
            arma::vec r1 = compute_residual(freq_table1_j, u1, weight.submat(0, j, n_sam_1 - 1, j));// weight.col(j).rows(0, n_sam_1 - 1)
            arma::vec r2 = compute_residual(freq_table2_j, u2, weight.submat(n_sam_1, j, n_sam - 1, j));// weight.col(j).rows(n_sam_1, n_sam - 1
            
            double sum_r1r2 = sum(r1 % r2); // sum of samples
            double sum_r1_squared_plus_r2_squared = sum(r1 % r1 + r2 % r2);
            
            // estiamte rho ; beta.n_rows
            double rho_est = ((1.0 / (n_sam - (K-1))) * sum_r1r2) / ((1.0 / (2 * n_sam - (K-1))) * sum_r1_squared_plus_r2_squared);// as j varies
            rho_estimates(j) = rho_est; // Store rho estimate for the jth taxon

            //DV = (freq_table.col(j) - u) % (weight.col(j)); //*: matrix multiplication; %: element-wise multiplication; freq_table.col
            DV1 = (freq_table1_j - u1) % weight.col(j).rows(0, n_sam_1 - 1) - rho_est * sqrt(v1 / v2) % weight.col(j).rows(n_sam_1, n_sam - 1) % (freq_table2_j - u2);
            DV2 = (freq_table2_j - u2) % weight.col(j).rows(n_sam_1, n_sam - 1) - rho_est * sqrt(v2 / v1) % weight.col(j).rows(0, n_sam_1 - 1) % (freq_table1_j - u1);

            S = X1.t() * DV1 + X2.t() * DV2;

            J_temp11 = - u1 % (1 - u1) % weight.col(j).rows(0, n_sam_1 - 1);
            J_temp21 = rho_est * sqrt(v1 / v2) % (u2 % (1 - u2)) % weight.col(j).rows(n_sam_1, n_sam - 1);
            J_temp22 = - u2 % (1 - u2) % weight.col(j).rows(n_sam_1, n_sam - 1);
            J_temp12 = rho_est * sqrt(v2 / v1) % (u1 % (1 - u1)) % weight.col(j).rows(0, n_sam_1 - 1);

           //J.eye(); //Use identity matrix to replace J.fill(0)
            J.fill(0);
            for (int i_sam = 0; i_sam < n_sam_1; i_sam++){
                 J = J + J_temp11(i_sam) * X1X1.slice(i_sam) + J_temp21(i_sam) * X2X1.slice(i_sam)
                       + J_temp22(i_sam) * X2X2.slice(i_sam) + J_temp12(i_sam) * X1X2.slice(i_sam);// J_temp11(i_sam) and similar variables are scalars
            }
            
            J_inv = inv(J);
            
            //std::cout << "success" << success << std::endl;

            //Firth
            if (prop_presence(j) < Firth_thresh){
                if (robust_var) {
                    Sigma.fill(0);
                    for (int i_sam = 0; i_sam < n_sam; i_sam++){
                        Sigma = Sigma + DV(i_sam) * DV(i_sam); //* XX.slice(i_sam);
                    }
                    Sigma = J_inv * Sigma * J_inv;
                } else {
                    Sigma = -J_inv;
                }
                J_temp11_1 = J_temp11 % (1 - 2*u1);
                J_temp21_1 = J_temp21 % (1 - 2*u2);
                J_temp22_1 = J_temp22 % (1 - 2*u2);
                J_temp12_1 = J_temp12 % (1 - 2*u1);

                for (int k = 0; k < K; k++) {

                  H_temp11 = J_temp11_1 % X.col(k).rows(0, n_sam_1 - 1);
                  H_temp21 = J_temp21_1 % X.col(k).rows(n_sam_1, n_sam - 1);
                  H_temp22 = J_temp22_1 % X.col(k).rows(n_sam_1, n_sam - 1);
                  H_temp12 = J_temp12_1 % X.col(k).rows(0, n_sam_1 - 1);

                    H.fill(0);
                    for (int i_sam = 0; i_sam < n_sam_1; i_sam++){
                        H = H + H_temp11(i_sam) * X1X1.slice(i_sam) + H_temp21(i_sam) * X2X1.slice(i_sam)
                              + H_temp22(i_sam) * X2X2.slice(i_sam) + H_temp12(i_sam) * X1X2.slice(i_sam);
                    }
                    S(k) = S(k) - 0.5*accu(H % Sigma); //trace(H * Sigma);
                }
            }
          
            // update
            //arma::vec lambda = {0.1, 0.1, 0.1};
            //arma::mat diag_lambda = diagmat(lambda);
            step =  0.1 * J_inv * S;
            
            if (sum(abs(step)) > 5*K){
                beta.col(j).fill(0);
            } else {
                beta.col(j) = beta.col(j) - step;
            }
            
            if (sum(abs(step) < tol) == K){
                break;
            }        
        } // iter
        final_step_sum(j) = sum(abs(step) < tol);

        if (get_var) {
            if (!robust_var) {
                Sigma.fill(0);
                for (int i_sam = 0; i_sam < n_sam; i_sam++){
                    Sigma = Sigma + DV(i_sam) * DV(i_sam); //* XX.slice(i_sam);
                }
                Sigma = J_inv * Sigma * J_inv;
            }
            variance.slice(j) = Sigma;
        }
    } // otu
    
    List L = List::create(Named("beta") = beta, Named("variance") = variance,
                          Named("rho_estimates") = rho_estimates, Named("final_step_sum") = final_step_sum);
    return L;
} // Newton()

// [[Rcpp::export]]
List  perm_Newton_cor(arma::mat freq_table, arma::mat Yr, 
                         arma::mat N_denominator, arma::mat X,
                         arma::mat beta_init, arma::mat weight,
                         arma::umat perm,
                         double tol, int iter_max, double Firth_thresh, bool robust_var,
                         arma::vec prop_presence, bool get_var) {
    int n_trait = Yr.n_cols;
    int n_otu = freq_table.n_cols;
    int n_perm = perm.n_cols; 
    arma::cube beta_est_temp(n_trait, n_otu, n_perm);
    arma::mat  final_step_converge(n_perm, n_otu);
    arma::mat  rho_estimates(n_perm, n_otu);
    arma::mat Yr_perm;
    arma::mat X_perm = X;
    uvec o; // a nonnegative integer
    
    for (int i_perm = 0; i_perm < n_perm; i_perm++){
        
        o = perm.col(i_perm);
        
        Yr_perm = Yr.rows(o);       
        X_perm.head_cols(n_trait) = Yr_perm;
        List res = Newton_cor(freq_table, N_denominator, X_perm, beta_init, weight, tol, iter_max, Firth_thresh, robust_var, prop_presence, get_var);
        arma::mat res_beta = res[0];
        arma::vec res_rho_estimates = res[2];
        arma::vec final_step_sum = res[3];
        beta_est_temp.slice(i_perm) = res_beta.head_rows(n_trait);
        rho_estimates.row(i_perm) = res_rho_estimates.t();
        final_step_converge.row(i_perm) = final_step_sum.t();
        
    }
    
    //return beta_est_temp;
     List L = List::create(Named("beta") = beta_est_temp,
                           Named("rho_estimates") = rho_estimates,
                           Named("final_step_converge") = final_step_converge);
     return L;
    
} // perm_Newton()






