#include "../20261118_experimental_counts_ar1_practical_01/ar1_kernel.cpp"

extern "C" int joint_set_histograms(void* handle, const double* histogram) {
  auto& d = *static_cast<Design*>(handle);
  for (int c = 0; c < 2; ++c) {
    std::copy_n(histogram + c*N, N, d.histogram[c].begin());
    if (std::abs(sum(d.histogram[c]) - 1.) > 1e-12) return 0;
  }
  d.cache_valid = false;
  return 1;
}

extern "C" int joint_live_states(void* handle, const double* p, int count,
                                 const int* pairs, double* values) {
  try {
    auto& d = *static_cast<Design*>(handle);
    if (!prepare(d, p, false)) return 0;
    for (int i = 0; i < count; ++i) {
      const auto& v = d.passages.at(pairs[2*i]).x.at(pairs[2*i+1]);
      if (!(sum(v)>0.)) return 0;
      std::copy(v.begin(), v.end(), values+i*N);
    }
    return 1;
  } catch (...) { return 0; }
}

extern "C" int joint_state_adjoint(void* handle, const double* p, int count,
                                   const int* pairs, const double* weights,
                                   double* grad) {
  try {
    auto& d = *static_cast<Design*>(handle);
    if (!prepare(d, p, false)) return 0;
    // The audited harvests are not used as donor ticks. The shared adjoint
    // clears donor ticks at entry and preserves the explicitly injected slots.
    for (int i = 0; i < count; ++i) {
      auto& seg = d.passages.at(pairs[2*i]);
      const int tick = pairs[2*i+1];
      if (seg.demand_index.at(tick) >= 0) return 0;
    }
    for (int i = 0; i < count; ++i) {
      auto& slot = d.passages[pairs[2*i]].external[pairs[2*i+1]];
      std::copy_n(weights+i*N, N, slot.begin());
    }
    const std::vector<double> zero(d.nslope);
    reverse(d, p, false, d.cached_rates, zero.data(), grad, false);
    for (int i = 0; i < count; ++i) {
      d.passages[pairs[2*i]].external[pairs[2*i+1]].fill(0.);
    }
    return 1;
  } catch (...) { return 0; }
}
