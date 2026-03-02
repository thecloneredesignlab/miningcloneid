#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_MATRIX(theta_hat);         // n x p local transformed parameters from per-sample fits
  DATA_VECTOR(sample_weight);     // length n
  DATA_SCALAR(tau_min);           // lower floor for hierarchical sd
  DATA_SCALAR(log_tau_prior_sd);  // weak prior scale on log_tau

  PARAMETER_VECTOR(mu);           // length p
  PARAMETER_VECTOR(log_tau);      // length p
  PARAMETER_MATRIX(theta_shrunk); // n x p MAP-shrunk local transformed parameters

  const int n = theta_hat.rows();
  const int p = theta_hat.cols();

  if (mu.size() != p) error("mu length must equal ncol(theta_hat)");
  if (log_tau.size() != p) error("log_tau length must equal ncol(theta_hat)");
  if (theta_shrunk.rows() != n || theta_shrunk.cols() != p) {
    error("theta_shrunk dimension must equal theta_hat dimension");
  }
  if (sample_weight.size() != n) error("sample_weight length must equal nrow(theta_hat)");

  Type nll = Type(0.0);
  vector<Type> tau(p);

  for (int j = 0; j < p; ++j) {
    tau(j) = tau_min + exp(log_tau(j));

    // Weak regularization on log_tau to stabilize estimates.
    nll += Type(0.5) * pow(log_tau(j) / log_tau_prior_sd, 2);

    for (int i = 0; i < n; ++i) {
      Type wi = sample_weight(i);
      if (!(wi > Type(0.0))) wi = Type(1.0);

      // Observation term: stay close to per-sample fitted theta_hat.
      Type d_obs = theta_hat(i, j) - theta_shrunk(i, j);
      nll += Type(0.5) * wi * d_obs * d_obs;

      // Hierarchical shrinkage toward shared center mu with scale tau.
      Type z = (theta_shrunk(i, j) - mu(j)) / tau(j);
      nll += Type(0.5) * z * z + log(tau(j));
    }
  }

  REPORT(mu);
  REPORT(tau);
  REPORT(theta_shrunk);
  ADREPORT(mu);
  ADREPORT(tau);

  return nll;
}
