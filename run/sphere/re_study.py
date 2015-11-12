from os import path
from numpy import exp, arange, pi, logspace, floor
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile


class test_case:
    def __init__(
        self, Re
    ):
        self.name = "re={0:0.2f}".format(Re)

        rho_f = 1.0                                                                                                           
        D = 1.0                                                                                                               
        R = D / 2                                                                                                             
        mu = 1.0 
        volume = pi / 6 * D**3
        V_t = Re * mu / rho_f / D
        g = 9.81
        rho_p = V_t * 18.0 * mu / D**2 / g + rho_f

        self.m = rho_p * volume
        self.ma = (rho_p - rho_f) * volume


def case_setup(ci):
    template_case = SolutionDirectory(
        "template", archive=None, paraviewLink=False)
    case = template_case.cloneCase("{0}".format(ci.name))

    stfproperties = ParsedParameterFile(path.join(
        case.name, "constant", "STFProperties"))

    stfproperties["apparentMass"] = ci.ma
    stfproperties["mass"] = ci.m
    stfproperties.writeFile()

Res = logspace(-2, 0, 10)
test_cases = [test_case(Re) for Re in Res]

if __name__ == "__main__":

    for tc in test_cases:
        case_setup(tc)
