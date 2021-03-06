from numpy import pi, logspace

"""
Based on Stokes drag law. It calculates and prints the terminal velocity a
Reynolds number based on it.

References:
  https://en.wikipedia.org/?title=Stokes'_law#Terminal_Velocity_of_Sphere_Falling_in_a_Fluid  

Notes:

Stokes' law applies when the Reynolds number, Re, of the particle is less than
0.1. Experimentally Stokes' law is found to hold within 1% for Re \leq 0.1,
within 3% for Re \leq 0.5 and within 9% Re \leq 1.0.[2] With increasing
Reynolds numbers, Stokes law begins to break down due to the increasing
importance of fluid inertia, requiring the use of empirical solutions to
calculate drag forces.
"""


rho_f = 1.0

D = 1.0
R = D / 2

mu = 1.0

g = 9.81

Res = logspace(-1, 0, 11)
volume = pi / 6 * D**3

for Re in Res:
    V_t = Re * mu / rho_f / D
    rho_p = V_t * 18.0 * mu / D**2 / g + rho_f
    # Need to account for apparent mass
    mass_p = rho_p * volume
    apparent_mass = (rho_p - rho_f) * volume
    print(
        'Re={0:g}; V_t={1:g} rhop={2:g} m={3:g} ma={4:g}, Fg={5:0.3g}'.format(
            Re, V_t, rho_p, mass_p, apparent_mass, apparent_mass * g))
