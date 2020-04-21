from sys import path
githubPath = r"\\arran.sms.ed.ac.uk\home\s1342398\GitHub"
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\editing')
path.append(githubPath + '\\interlaced_model_creation\\editing\\3D\\Mechanical Properties\\Parametric')
import interlacedRVEImplicit as intRVE
from itertools import product
from abaqusConstants import *
import csv


def parametricStudy():
    # define parameters and create design space
    widths = (15.0, 20.0, 25.0)
    spacing = (1, 2, 3, 4)
    designSpace = product(widths, spacing)

    # define non parameter variables
    tapeAngles = (0, 90)
    tapeThickness = 0.18
    undulationRatio = 0.18
    displacements = [0.1, UNSET, UNSET]

    reportFilePath = 'C:\\Workspace\\3D_RVE\\Parametric\\PS_Summary.csv'
    headings = ['tapeWidth', 'tapeSpacing', 'E11']
    stiffnessData = {}

    # run simulations
    with open(reportFilePath, 'ab') as csvfile:
        csvWriter = csv.writer(csvfile, delimiter=',',
                               quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvWriter.writerow(headings)                       
        for study in designSpace:
            try:
                tapeWidths = (study[0], )* len(tapeAngles)
                tapeSpacing = study[1]
                studyStiffness = intRVE.main(tapeAngles, tapeWidths,
                                             tapeSpacing, tapeThickness,
                                             undulationRatio, displacements)
                stiffnessData[study] = studyStiffness
                csvWriter.writerow([tapeWidths[0], tapeSpacing, studyStiffness])
            except:
                csvWriter.writerow([tapeWidths[0], 'NaN', 'NaN'])
                continue


if __name__ == '__main__':
    parametricStudy()