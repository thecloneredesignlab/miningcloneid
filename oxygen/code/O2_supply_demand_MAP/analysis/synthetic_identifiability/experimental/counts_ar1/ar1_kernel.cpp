#include "../20261118_experimental_counts_practical_01/count_kernel.cpp"

extern "C" int ar1_count_loss(void* handle, const double* p, int n, const int* slots,
                              const int* previous, const int* gap, const double* target,
                              double* out, double* grad) {
  try {
    auto& d = *static_cast<Design*>(handle);
    if (!prepare(d, p, false)) return 0;
    const double rho = p[NC];
    std::vector<double> residual(n), weighted(n), source(d.nslope);
    double q = 0., logdet = 0., dq = 0., dlogdet = 0.;
    for (int i = 0; i < n; ++i) residual[i] = d.logc[slots[i]] - target[i];
    for (int i = 0; i < n; ++i) {
      const int j = previous[i];
      if (j < 0) {
        q += residual[i] * residual[i];
        weighted[i] += residual[i];
      } else {
        const double a = std::pow(rho, gap[i]);
        const double da = gap[i] == 1 ? 1. : gap[i] * std::pow(rho, gap[i]-1);
        const double v = 1.-a*a;
        const double e = residual[i] - a*residual[j];
        const double z = e/v;
        q += e*z;
        logdet += std::log1p(-a*a);
        dq += 2.*z*da*(-residual[j]+a*z);
        dlogdet -= 2.*a*da/v;
        weighted[i] += z;
        weighted[j] -= a*z;
      }
    }
    const double factor = std::exp(logdet/n);
    out[0] = .5*q*factor;
    out[1] = *std::min_element(d.margins.begin(),d.margins.end());
    out[2] = q;
    out[3] = logdet;
    if (grad) {
      for (int i = 0; i < n; ++i) {
        const int j = slots[i];
        if (d.cells[j]>1.) source[j] = factor * weighted[i]*scale/d.cells[j];
      }
      reverse(d,p,false,d.cached_rates,source.data(),grad,false);
      grad[NC] = .5*factor*(dq+q*dlogdet/n);
    }
    return std::isfinite(out[0]);
  } catch (...) {
    return 0;
  }
}
