
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

# include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ComputeAlphaBeta(arma::vec y, arma::vec x, arma::mat WW,
                      arma::vec weight, arma::mat Z, arma::vec G) {

    // Inputs
    // y, x: n x 1
    // WW: n x k (controls + intercept)
    // weight: n x 1 (non-negative weights)
    // Z: n x m (local shares/instruments)
    // G: m x 1 (global shocks)

    const int n = x.n_elem;

    // Guard: if weight length mismatches, fall back to uniform weights
    if (static_cast<int>(weight.n_elem) != n) {
        weight.set_size(n);
        weight.ones();
    }

    // Precompute weighted design cross-products: Q = W' * diag(w) * W
    // Achieved efficiently via column-wise scaling by w
    arma::mat WW_w = WW.each_col() % weight;           // n x k
    arma::mat Q = WW.t() * WW_w;                       // k x k

    // Solve for projection coefficients without forming the residual-maker
    // For a vector v: fitted_v = WW * solve(Q, WW' * (w % v))
    auto resid_vec = [&](const arma::vec& v) -> arma::vec {
        arma::vec rhs = WW.t() * (weight % v);         // k x 1
        arma::vec coef = arma::solve(Q, rhs, arma::solve_opts::fast);
        return v - WW * coef;                          // n x 1 residuals
    };

    arma::vec xx = resid_vec(x);
    arma::vec yy = resid_vec(y);

    // Compute cross-products efficiently
    arma::vec wx = weight % xx;                        // n x 1
    arma::vec wy = weight % yy;                        // n x 1
    arma::vec zt_w_x = Z.t() * wx;                     // m x 1
    arma::vec zt_w_y = Z.t() * wy;                     // m x 1

    // Rotemberg weights (alpha)
    arma::vec num = G % zt_w_x;                        // m x 1
    double denom = arma::as_scalar(G.t() * zt_w_x);    // scalar
    arma::vec Alpha = num / denom;                     // m x 1

    // Just-identified coefficients (beta)
    arma::vec Beta = zt_w_y / zt_w_x;                  // m x 1 (element-wise)

    return List::create(Alpha, Beta);
}
