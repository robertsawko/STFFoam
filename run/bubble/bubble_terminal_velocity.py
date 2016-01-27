from numpy import sqrt, linspace
from scipy.optimize import brentq

"""
The calculations are heavily based on:

    Tomiyama, A, Kataoka, I., Zun, I. and Sakaguchi, T. "Drag
    Coefficients of single bubbles under normal and micro gravity
    conditions", JSME International Journal, Series B Vol 41, No 2 1998
"""


def EotvosNumber(rhoc, rhod, sigma, d, geff=9.81):
    """Eotvos number expression

    Non-dimensional number that measures the importance of surface
    tension compared to body forces.
    """
    return geff * abs(rhod - rhoc) * d**2 / sigma


def MortonNumber(rhoc, rhod, sigma, muc, geff=9.81):
    """Morton number expression

    Non-dimensional number chracterising fluid material properties
    """
    return \
        geff * muc**4 * abs(rhod - rhoc) /\
        (sigma**3 * rhoc**2)


def Cd(Re, Eo):
    """Tomiyama drag coefficient for pure system
    """
    return max(
        min(
            16 / Re * (1 + 0.15 * Re**0.687),
            48 / Re),
        8 / 3 * Eo / (Eo + 4))


def ReEoMo_equation(Re, Eo, Mo):
    """Force balance equation in nondimensional form

    Tomiyama, A, Kataoka, I., Zun, I. and Sakaguchi, T. "Drag
    Coefficients of single bubbles under normal and micro gravity
    conditions", JSME International Journal, Series B Vol 41, No 2 1998
    """
    return Re**2 * Cd(Re, Eo) - 4 / 3 * sqrt(Eo**3 / Mo)


def terminal_velocity(Re, muc, rhoc, d):
    return Re * muc / d / rhoc

rhow = 1000
rhob = 1.48e-05

D = 1.0
R = D / 2

nuw = 1.0e-6
muw = rhow * nuw

sigma = 0.0707106
diameters = linspace(0.001, 0.005, 5)

for d in diameters:
    Eo = EotvosNumber(rhoc=rhow, rhod=rhob, sigma=sigma, d=d)
    Mo = MortonNumber(rhoc=rhow, rhod=rhob, sigma=sigma, muc=muw)
    result = brentq(lambda x: ReEoMo_equation(x, Eo, Mo), 0.001, 2000)
    print('{0:0.3f} {1:0.3f}'.format(
        d, terminal_velocity(result, muw, rhow, d)))
