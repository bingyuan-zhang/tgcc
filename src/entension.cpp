#include <RcppArmadillo.h>
#include "extension.h"

// [[Rcpp::depends(RcppArmadillo)]]

double loss_func_fista(const arma::colvec & y,
                       const arma::colvec & x,
                       double lam) {
  double term1 = pow(arma::norm(y - x, 2), 2) / 2 ;
  double term2 = lam * arma::norm(x, 2);
  return term1 + term2;
}

arma::colvec fista_descent(const arma::colvec & y,
  double lam,
  int maxiter) {
  int n = y.n_elem;
  arma::colvec x_old(n);
  arma::colvec x_new(n);
  arma::colvec v_old(n);
  arma::colvec v_new(n);
  arma::colvec u_tmp(n);
  arma::colvec uu(n);
  arma::colvec u(n);
  arma::colvec loss(maxiter);

  double loss_old = 1e5;
  double t;
  double thresh;
  double loss_new;
  bool descent;
  int no_descent_tol = 10;
  int no_descent = 0;

  for (int iter = 0; iter < maxiter; ++iter) {

    // FISTA step
    t = 2.0 / (iter + 2);

    uu = (1 - t) * x_old + t * v_old;
    u_tmp = uu - 0.5 * (uu - y);
    thresh = std::max(0.0, 1.0 - lam / arma::norm(u_tmp,2) );
    u = u_tmp * thresh;

    loss_new = loss_func_fista(y, u, lam);

    descent = loss_new < loss_old;
    if (descent) {
      x_new = u;
      loss_old = loss_new;
      no_descent = 0;
    } else {
      x_new = x_old;
      no_descent++;
    }

    v_new = x_old + 1 / t * (u - x_old);

    if (no_descent >= no_descent_tol) break;

    x_old = x_new;
    v_old = v_new;
    loss[iter] = loss_old;
  }

  return x_old;
}

double calculate_sparseclustering_loss(const arma::mat& data,
                                       const arma::mat& U,
                                       double lam,
                                       double gam) {

  double term1 = 0.5 * arma::norm(data - U, "fro") * arma::norm(data - U, "fro");
  int n = data.n_rows;
  int p = data.n_cols;
  double term2 = 0;
  double term3 = 0;

  // term 2
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      term2 += sum(abs(U.col(i) - U.col(j)));
    }
  }
  // term 3
  for (int i = 0; i < n; ++i) {
    term3 += norm(U.row(i), 2);
  }

  return (term1 + lam*term2 + gam*term3);
}

// [[Rcpp::export]]
arma::mat updateUSC(const arma::mat& U,
                    const arma::mat& U0,
                    double lam,
                    double gam,
                    int max_iter,
                    double precision_thres,
                    const std::vector<double>& V,
                    const std::vector<double>& Type,
                    const std::vector<double>& Pa,
                    const std::vector<double>& WV,
                    const std::vector<double>& WE,
                    const std::vector<std::vector<double>>& C) {

  int n = U.n_rows;
  int p = U.n_cols;

  arma::mat Unew = U0;
  arma::mat P = arma::zeros(n, p);
  arma::mat Q = arma::zeros(n, p);
  arma::mat Y = arma::zeros(n, p);

  arma::vec diff_iter = arma::zeros(max_iter);
  arma::vec loss_iter = arma::zeros(max_iter);

  double diff = arma::datum::inf;
  double loss_old = calculate_sparseclustering_loss(U0, U0, lam, gam);
  double loss_cur;

  std::vector<double> tempcol(n);
  std::vector<double> tempcolres(n);

  for (int iter = 0; iter < max_iter; ++iter) {

    // Update Y
    for (int j = 0; j < p; ++j) {
      Y.col(j) = fista_descent(Unew.col(j) + P.col(j), gam, max_iter);

      P.col(j) = Unew.col(j) + P.col(j) - Y.col(j);

      tempcol = arma::conv_to < std::vector<double> > ::from(Y.col(j) + Q.col(j));
      tempcolres = dp (tempcol, lam, V, Type, Pa, WV, WE, C);
      Unew.col(j) = arma::conv_to< arma::colvec > ::from(tempcolres);

      Q.col(j) = Y.col(j) + Q.col(j) - Unew.col(j);
    }

    // Calculate current loss and check convergence
    loss_cur = calculate_sparseclustering_loss(U0, Unew, lam, gam);

    diff = std::abs(loss_cur - loss_old) / loss_old;
    loss_iter(iter) = loss_cur;
    diff_iter(iter) = diff;
    loss_old = loss_cur;

    if (diff < precision_thres) {
      break;
    }

  }

  return Unew;
}


double calculate_biclustering_loss(const arma::mat& data,
                                   const arma::mat& U,
                                   double lam,
                                   double gam) {

  double term1 = 0.5 * arma::norm(data - U, "fro") * arma::norm(data - U, "fro");
  int n = data.n_rows;
  int p = data.n_cols;
  double term2 = 0;
  double term3 = 0;

  // term 2
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      term2 += sum(abs(U.col(i) - U.col(j)));
    }
  }
  // term 3
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      term3 += sum(abs(U.row(i) - U.row(j)));
    }
  }

  return (term1 + lam*term2 + gam*term3);
}

// [[Rcpp::export]]
arma::mat updateUBC(const arma::mat& U, const arma::mat& U0,
                    double lam, double gam, int max_iter, double precision_thres,
                    const std::vector<double>& V_s, const std::vector<double>& Type_s, const std::vector<double>& Pa_s,
                    const std::vector<double>& WV_s, const std::vector<double>& WE_s, const std::vector<std::vector<double>>& C_s,
                    const std::vector<double>& V_f, const std::vector<double>& Type_f, const std::vector<double>& Pa_f,
                    const std::vector<double>& WV_f, const std::vector<double>& WE_f, const std::vector<std::vector<double>>& C_f) {

  int n = U.n_rows;
  int p = U.n_cols;

  arma::mat Unew = U;
  arma::mat P = arma::zeros(n, p);
  arma::mat Q = arma::zeros(n, p);
  arma::mat Y = arma::zeros(n, p);

  arma::vec diff_iter = arma::zeros(max_iter);
  arma::vec loss_iter = arma::zeros(max_iter);

  double diff = arma::datum::inf;
  double loss_old = calculate_biclustering_loss(U0, U0, lam, gam);
  double loss_cur;
  std::vector<double> tempcol(n);
  std::vector<double> tempcolres(n);
  std::vector<double> temprow(p);
  std::vector<double> temprowres(p);


  for (int iter = 0; iter < max_iter; ++iter) {

    // Update Y
    for (int j = 0; j < p; ++j) {
      tempcol = arma::conv_to < std::vector<double> > ::from(Unew.col(j) + P.col(j));
      tempcolres = dp (tempcol, lam, V_s, Type_s, Pa_s, WV_s, WE_s, C_s);
      Y.col(j) = arma::conv_to< arma::colvec >::from(tempcolres);
    }
    // Update P
    P = Unew + P - Y;

    // Update U
    for (int j = 0; j < n; ++j) {
      temprow = arma::conv_to<std::vector<double>> ::from(Y.row(j) + Q.row(j));
      temprowres = dp (temprow, gam, V_f, Type_f, Pa_f, WV_f, WE_f, C_f);
      Unew.row(j) = arma::conv_to< arma::rowvec >::from(temprowres);
    }

    // Update Q
    Q = Y + Q - Unew;

    // Calculate current loss and check convergence
    loss_cur = calculate_biclustering_loss(U0, Unew, lam, gam);
    diff = std::abs((loss_cur - loss_old) / loss_old);

    loss_iter(iter) = loss_cur;
    diff_iter(iter) = diff;
    loss_old = loss_cur;

    if (diff < precision_thres) break;
  }
  return Unew;
}
