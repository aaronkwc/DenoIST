#include <Rcpp.h>

using namespace Rcpp;

// Poisson mixture EM core loop.
// [[Rcpp::export]]
List poisson_mixture_em_cpp(NumericVector x,
                            NumericVector s,
                            double lambda1,
                            double lambda2,
                            double pi,
                            int max_iter,
                            double tol) {
  const R_xlen_t n = x.size();
  NumericVector gamma(n);

  auto log_likelihood = [&](double l1, double l2, double p) {
    double acc = 0.0;
    for (R_xlen_t i = 0; i < n; ++i) {
      const double t1 = p * R::dpois(x[i], s[i] * l1, false);
      const double t2 = (1.0 - p) * R::dpois(x[i], s[i] * l2, false);
      const double mix = t1 + t2;
      acc += std::log(mix);
    }
    return acc;
  };

  double log_lik = log_likelihood(lambda1, lambda2, pi);

  for (int iter = 0; iter < max_iter; ++iter) {
    double gamma_s = 0.0;
    double one_minus_gamma_s = 0.0;
    double gamma_x = 0.0;
    double one_minus_gamma_x = 0.0;
    double gamma_sum = 0.0;

    for (R_xlen_t i = 0; i < n; ++i) {
      const double t1 = pi * R::dpois(x[i], s[i] * lambda1, false);
      const double t2 = (1.0 - pi) * R::dpois(x[i], s[i] * lambda2, false);
      const double denom = t1 + t2;
      double g = (denom > 0.0) ? (t1 / denom) : 0.5;
      if (!R_finite(g)) {
        g = 0.5;
      }
      gamma[i] = g;
      gamma_s += g * s[i];
      one_minus_gamma_s += (1.0 - g) * s[i];
      gamma_x += g * x[i];
      one_minus_gamma_x += (1.0 - g) * x[i];
      gamma_sum += g;
    }

    if (gamma_s <= 0.0 || one_minus_gamma_s <= 0.0) {
      break;
    }

    lambda1 = gamma_x / gamma_s;
    lambda2 = one_minus_gamma_x / one_minus_gamma_s;
    pi = gamma_sum / static_cast<double>(n);

    const double new_log_lik = log_likelihood(lambda1, lambda2, pi);
    if (!R_finite(new_log_lik) || std::abs(new_log_lik - log_lik) < tol) {
      break;
    }
    log_lik = new_log_lik;
  }

  return List::create(_["lambda1"] = lambda1,
                      _["lambda2"] = lambda2,
                      _["pi"] = pi,
                      _["log_lik"] = log_lik,
                      _["gamma"] = gamma);
}

// Mini-batch Poisson mixture EM core loop.
// [[Rcpp::export]]
List poisson_mixture_batch_em_cpp(NumericVector x,
                                  NumericVector s,
                                  IntegerVector non_zero_indices,
                                  int batch_size,
                                  int n_epochs,
                                  double tau,
                                  double alpha,
                                  double lambda1,
                                  double lambda2,
                                  double pi) {
  const R_xlen_t n_filtered = non_zero_indices.size();

  for (int epoch = 0; epoch < n_epochs; ++epoch) {
    double acc_gamma_x = 0.0;
    double acc_gamma_s = 0.0;
    double acc_gamma = 0.0;
    double acc_1mg_x = 0.0;
    double acc_1mg_s = 0.0;

    for (R_xlen_t start = 0; start < n_filtered; start += batch_size) {
      const R_xlen_t end = std::min(start + static_cast<R_xlen_t>(batch_size), n_filtered);
      for (R_xlen_t i = start; i < end; ++i) {
        const int idx = non_zero_indices[i] - 1;
        const double xi = x[idx];
        const double si = s[idx];
        const double t1 = pi * R::dpois(xi, si * lambda1, false);
        const double t2 = (1.0 - pi) * R::dpois(xi, si * lambda2, false);
        const double denom = t1 + t2;
        double g = (denom > 0.0) ? (t1 / denom) : 0.5;
        if (!R_finite(g)) {
          g = 0.5;
        }
        acc_gamma_x += g * xi;
        acc_gamma_s += g * si;
        acc_gamma += g;
        acc_1mg_x += (1.0 - g) * xi;
        acc_1mg_s += (1.0 - g) * si;
      }
    }

    if (acc_gamma_s <= 0.0 || acc_1mg_s <= 0.0) {
      break;
    }

    const double lambda1_new = acc_gamma_x / acc_gamma_s;
    const double lambda2_new = acc_1mg_x / acc_1mg_s;
    const double pi_new = acc_gamma / static_cast<double>(n_filtered);

    const double learning_rate = (tau > 0.0) ? (1.0 / (1.0 + tau / static_cast<double>(n_filtered))) : 1.0;

    if (alpha > 0.0 && epoch > 0) {
      lambda1 = alpha * lambda1 + (1.0 - alpha) * learning_rate * lambda1_new + (1.0 - learning_rate) * lambda1;
      lambda2 = alpha * lambda2 + (1.0 - alpha) * learning_rate * lambda2_new + (1.0 - learning_rate) * lambda2;
      pi = alpha * pi + (1.0 - alpha) * learning_rate * pi_new + (1.0 - learning_rate) * pi;
    } else {
      lambda1 = learning_rate * lambda1_new + (1.0 - learning_rate) * lambda1;
      lambda2 = learning_rate * lambda2_new + (1.0 - learning_rate) * lambda2;
      pi = learning_rate * pi_new + (1.0 - learning_rate) * pi;
    }
  }

  double log_lik = 0.0;
  for (R_xlen_t i = 0; i < n_filtered; ++i) {
    const int idx = non_zero_indices[i] - 1;
    const double t1 = pi * R::dpois(x[idx], s[idx] * lambda1, false);
    const double t2 = (1.0 - pi) * R::dpois(x[idx], s[idx] * lambda2, false);
    const double mix = t1 + t2;
    if (mix <= 0.0) {
      log_lik = R_NegInf;
      break;
    }
    log_lik += std::log(mix);
  }

  return List::create(_["lambda1"] = lambda1,
                      _["lambda2"] = lambda2,
                      _["pi"] = pi,
                      _["log_lik"] = log_lik);
}
