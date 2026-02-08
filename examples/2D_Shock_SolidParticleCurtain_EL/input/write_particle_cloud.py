import numpy as np

# -----------------------------
# User parameters
# -----------------------------
N = 100                    # number of particles
x_min, x_max = 0.01, 0.02    # curtain thickness in x (m)
y_min, y_max = -0.125, 0.125   # height in y (m)
z_val = 0.0                  # 2D simulation

r_mean = 10e-6               # mean particle radius (m)
r_spread = 0.2               # ±20% polydispersity

output_file = "lag_particles.dat"

# -----------------------------
# Generate particle properties
# -----------------------------
rng = np.random.default_rng(seed=1)  # reproducible

x = rng.uniform(x_min, x_max, N)
y = rng.uniform(y_min, y_max, N)
z = np.full(N, z_val)

u = np.zeros(N)
v = np.zeros(N)
w = np.zeros(N)

radius = r_mean * (1.0 + r_spread * (rng.random(N) - 0.5))
interface_velocity = np.zeros(N)

# -----------------------------
# Stack and write file
# -----------------------------
data = np.column_stack([
    x, y, z,
    u, v, w,
    radius,
    interface_velocity
])

np.savetxt(
    output_file,
    data,
    fmt="%.8e",
    comments=""
)

print(f"Wrote {N} particles to '{output_file}'")
