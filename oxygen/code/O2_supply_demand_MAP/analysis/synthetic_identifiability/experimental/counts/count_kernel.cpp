#include "../20261118_experimental_OLS_practical_01/ols_kernel.cpp"

extern "C" int count_loss(void* handle, const double* p, int gaussian, int n,
                          const int* slots, const double* target, const double* precision,
                          int dense, double* out, double* grad) {
  try {
    auto& d = *static_cast<Design*>(handle);
    if (!prepare(d, p, gaussian != 0)) return 0;
    std::vector<double> residual(n), weighted(n), source(d.nslope);
    for (int j = 0; j < n; ++j) residual[j] = d.logc[slots[j]] - target[j];
    out[0] = 0.;
    for (int j = 0; j < n; ++j) {
      if (dense) {
        for (int k = 0; k < n; ++k) weighted[j] += precision[j*n+k]*residual[k];
      } else {
        weighted[j] = precision[j]*residual[j];
      }
      out[0] += .5*weighted[j]*residual[j];
      const int slot = slots[j];
      if (grad && d.cells[slot] > 1.) source[slot] = scale/d.cells[slot]*weighted[j];
    }
    out[1] = *std::min_element(d.margins.begin(), d.margins.end());
    if (grad) reverse(d, p, gaussian != 0, d.cached_rates, source.data(), grad, false);
    return 1;
  } catch (...) {
    return 0;
  }
}
