from os import path
from numpy import linspace, floor
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParameterFile import ParameterFile


class test_case:
    def __init__(self, d):
        self.bubble_diameter = d
        self.width = 8 * d
        self.height = 10 * d
        self.n = 50

        self.name = "{0:0.0f}mm_bubble".format(d*1000)


def case_setup(ci):
    template_case = SolutionDirectory(
        "template", archive=None, paraviewLink=False)
    case = template_case.cloneCase("{0}".format(ci.name))

    blockmeshcfg = ParameterFile(path.join(
        case.name, "system", "blockMeshDict.cfg"))

    blockmeshcfg.replaceParameter(
        "diameter", '{0:f}'.format(ci.bubble_diameter))
    blockmeshcfg.replaceParameter("height", '{0:f}'.format(ci.height))
    blockmeshcfg.replaceParameter("width", '{0:f}'.format(ci.width))
    blockmeshcfg.replaceParameter("nside", '{0:d}'.format(ci.n))
    blockmeshcfg.writeFile()


diameters = linspace(0.001, 0.005, 5)
test_cases = [test_case(d) for d in diameters]

if __name__ == "__main__":
    for tc in test_cases:
        case_setup(tc)
