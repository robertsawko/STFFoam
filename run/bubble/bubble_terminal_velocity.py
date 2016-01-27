from numpy import sqrt, linspace
from scipy.optimize import brentq

"""
Based on Tomiyama
"""

geff = 9.81


def Eof(rhoc, rhod, sigma, d):
    return geff * abs(rhod - rhoc) * d**2 / sigma


def Mof(rhoc, rhod, sigma, muc):
    return \
        geff * muc**4 * abs(rhod - rhoc) /\
        (sigma**3 * rhoc**2)


def Cd(Re, Eo):
    return max(
        min(
            16 / Re * (1 + 0.15 * Re**0.687),
            48 / Re),
        8 / 3 * Eo / (Eo + 4))


def ReEoMo_equation(Re, Eo, Mo):
    return Re**2 * Cd(Re, Eo) - 4 / 3 * sqrt(Eo**3 / Mo)


def terminal_velocity(Re, muc, rhoc, d):
    return Re * muc / d / rhoc

rhow = 1000
rhob = 1.48e-05

D = 1.0
R = D / 2

nuw = 1.0e-6
nub = 1.48e-5
muw = rhow * nuw

sigma = 0.0707106
ds = linspace(0.001, 0.005, 5)

for d in ds:
    Eo = Eof(rhoc=rhow, rhod=rhob, sigma=sigma, d=d)
    Mo = Mof(rhoc=rhow, rhod=rhob, sigma=sigma, muc=muw)
    result = brentq(lambda x: ReEoMo_equation(x, Eo, Mo), 0.001, 2000)
    print('{0:0.3f} {1:0.3f}'.format(
        d, terminal_velocity(result, muw, rhow, d)))
