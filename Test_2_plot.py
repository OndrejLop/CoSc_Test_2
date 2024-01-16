import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define the damped harmonic oscillator equation
def damped_harmonic_oscillator(y, t, m, c, k):
    x, v = y
    dxdt = v
    dvdt = -(c/m) * v - (k/m) * x
    return [dxdt, dvdt]

# Set initial conditions and parameters
y0 = [1.0, 0.0]  # Initial displacement and velocity
m = 1.0          # Mass
c = 0.5          # Damping coefficient
k = 10.0         # Spring constant

# Set time points for the solution
t = np.linspace(0, 10, 1000)

# Solve the differential equation
solution = odeint(damped_harmonic_oscillator, y0, t, args=(m, c, k))

# Extract displacement and velocity from the solution
x, v = solution.T

# Plot the results
plt.plot(t, x, label='Displacement')
plt.plot(t, v, label='Velocity')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

# Save the plot as a PNG file
plt.savefig('damped_harmonic_oscillator_plot.png')

# Display the plot
plt.show()
