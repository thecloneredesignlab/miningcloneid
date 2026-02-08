import time
# Import your actual function
from PujanEarlyVersionModel import ploidy_forcast

# Dummy data
test_ploidy = {2.0: 1e9, 3.0: 0, 4.0: 0.5e9}
drug = "gemcitabine"

# Time one simulation
t0 = time.time()
# Run the exact function your loop calls
ploidy_forcast(test_ploidy, drug, 7.0, N_SIMS=1000, R_BASE=0.575)
t1 = time.time()

elapsed = t1 - t0
print(f"One simulation took: {elapsed:.4f} seconds")

# Calculate prediction for Depth 10 (on 5 cores)
total_hours = (9765625 * elapsed) / 5 / 3600
print(f"Estimated runtime for Depth 10: {total_hours:.2f} hours")