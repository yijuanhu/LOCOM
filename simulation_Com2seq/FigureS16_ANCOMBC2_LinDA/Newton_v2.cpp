#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube CalculateXX(const arma::mat & X){
    
    int K = X.n_cols;
    int n_sam = X.n_rows;
    arma::cube XX(K, K, n_sam);
    
    for (int s = 0; s < n_sam; s++){
        for (int i = 0; i < K; i++){
            XX(i, i, s) = X(s, i) * X(s, i);
            for (int j = i+1; j < K; j++){
                XX(i, j, s) = X(s, i) * X(s, j);
                XX(j, i, s) = XX(i, j, s);
            }
        }
    }
    return XX;
}


arma::cube UpdateXX(const arma::cube & XX, const arma::mat & X, arma::mat Yr){
    
    int n_sam = XX.n_slices;
    int n_trait = Yr.n_cols;
    int K = XX.n_rows;
    arma::cube XX_perm = XX;
    
    for (int s = 0; s < n_sam; s++){
        XX_perm.slice(s).submat(0, 0, n_trait-1, n_trait-1) = trans(Yr.row(s)) * Yr.row(s);
        XX_perm.slice(s).submat(0, n_trait, n_trait-1, K-1) = trans(Yr.row(s)) * X.submat(s, n_trait, s, K-1);
        XX_perm.slice(s).submat(n_trait, 0, K-1, n_trait-1) = trans(XX_perm.slice(s).submat(0, n_trait, n_trait-1, K-1));
    }
    return XX_perm;
}

// [[Rcpp::export]]
List Newton(arma::mat freq_table, arma::mat X, arma::cube XX,
           arma::mat beta_init, arma::mat weight,
           double tol, int iter_max, double Firth_thresh, bool robust_var,
           arma::vec prop_presence, bool get_var) {

    int n_otu = freq_table.n_cols;
    int n_sam = freq_table.n_rows;
    int K = X.n_cols;

    arma::mat beta = beta_init;
    arma::vec z;
    arma::vec u;
    arma::vec DV; 
    arma::vec S;
    arma::vec J_temp;
    arma::vec J_temp_1;
    arma::mat J(K, K, fill::zeros);
    arma::mat J_inv;
    arma::mat Sigma(K, K, fill::zeros);
    arma::vec H_temp;
    arma::mat H(K, K, fill::zeros);
    arma::vec step;
    arma::cube variance;
    if (get_var) variance = cube(K, K, n_otu, fill::zeros);
    arma::vec final_step_sum(n_otu);
    
    for (int j = 0; j < n_otu; j++){

        beta.col(j) = beta_init.col(j);

        for (int i = 0; i < iter_max; i++){

            z = exp(X * beta.col(j));
            u = z / (1 + z);
            DV = (freq_table.col(j) - u) % (weight.col(j)); //*: matrix multiplication; %: element-wise multiplication
            S = X.t() * DV;
            
            J_temp = - (u % (1 - u)) % (weight.col(j));

           //J.eye(); //Use identity matrix to replace J.fill(0)
            J.fill(0);
            for (int i_sam = 0; i_sam < n_sam; i_sam++){
                J = J + J_temp(i_sam) * XX.slice(i_sam);
            }
            J_inv = inv(J);

            //Firth

            if (prop_presence(j) < Firth_thresh){
                
                if (robust_var) {
                    Sigma.fill(0);
                    for (int i_sam = 0; i_sam < n_sam; i_sam++){
                        Sigma = Sigma + DV(i_sam) * DV(i_sam) * XX.slice(i_sam);
                    }
                    Sigma = J_inv * Sigma * J_inv;
                } else {
                    Sigma = -J_inv;
                }
                
                J_temp_1 = J_temp % (1 - 2*u);

                for (int k = 0; k < K; k++) {

                    H_temp = J_temp_1 % X.col(k);

                    H.fill(0);
                    for (int i_sam = 0; i_sam < n_sam; i_sam++){
                        H = H + H_temp(i_sam) * XX.slice(i_sam);
                    }
                    S(k) = S(k) - 0.5*accu(H % Sigma); //trace(H * Sigma);
                }
            }

            // update
            step =  J_inv * S;

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
                    Sigma = Sigma + DV(i_sam) * DV(i_sam) * XX.slice(i_sam);
                }
                Sigma = J_inv * Sigma * J_inv;
            }
            variance.slice(j) = Sigma;
        }
        
        // Rcout << "otu: " << j << " n_iter:" << i << "\n";
    } // otu
    
    List L = List::create(Named("beta") = beta, Named("variance") = variance,
                          Named("final_step_sum") = final_step_sum);

    return L;
} // Newton()

// [[Rcpp::export]]
List perm_Newton(arma::mat freq_table, arma::mat Yr, arma::mat X, arma::cube XX,
                         arma::mat beta_init, arma::mat weight,
                         arma::umat perm,
                         double tol, int iter_max, double Firth_thresh, bool robust_var,
                         arma::vec prop_presence, bool get_var) {
    int n_trait = Yr.n_cols;
    int n_otu = freq_table.n_cols;
    int n_perm = perm.n_cols; 
    
    arma::cube beta_est_temp(n_trait, n_otu, n_perm);
    arma::mat  final_step_converge(n_perm, n_otu);

    arma::mat Yr_perm;
    arma::mat X_perm = X;
    uvec o;
    
    
    for (int i_perm = 0; i_perm < n_perm; i_perm++){
        
        o = perm.col(i_perm);
        
        Yr_perm = Yr.rows(o);             
        X_perm.head_cols(n_trait) = Yr_perm;
        XX = UpdateXX(XX, X, Yr_perm);
        
        List res = Newton(freq_table, X_perm, XX, beta_init, weight, tol, iter_max, Firth_thresh, robust_var, prop_presence, get_var);
        arma::mat res_beta = res[0];
        arma::vec final_step_sum = res[2];
        beta_est_temp.slice(i_perm) = res_beta.head_rows(n_trait);
        final_step_converge.row(i_perm) = final_step_sum.t();
    }
    
    //return beta_est_temp;
    
    List L = List::create(Named("beta") = beta_est_temp,
                          Named("final_step_converge") = final_step_converge);

    return L;
    
} // perm_Newton()






