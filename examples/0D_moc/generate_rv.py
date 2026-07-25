import numpy as np 
import numpy as np

# Reproducible random number generator
rng = np.random.default_rng(seed=42)

# Distribution parameters
mean = np.array([1.0, 0.0])
cov = np.array([
    [0.01, 0.0],
    [0.0, 0.01]
])

# Generate 1000 samples
samples = rng.multivariate_normal(mean, cov, size=1000)

# Split into radius and velocity samples
rad = samples[:, 0]
vel = samples[:, 1]

# Compute sample statistics
sample_mean = samples.mean(axis=0)
sample_cov = np.cov(samples, rowvar=False)

print("Sample mean:")
print(sample_mean)

print("\nSample covariance:")
print(sample_cov)

# Save to text files
np.savetxt("Rad.txt", rad, fmt="%.16e")
np.savetxt("Vel.txt", vel, fmt="%.16e")

print("\nSaved:")
print("  Rad.txt")
print("  Vel.txt")