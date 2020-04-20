from sys import path
githubPath = 'C:\\Users\\rutge\\Documents\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\editing')
path.append(githubPath + '\\interlaced_model_creation\\editing\\3D\\Mechanical Properties\\Parametric')
import interlacedRVEImplicit
from itertools import product
from abaqusConstants import *
import csv


def parametricStudy():
    # define parameters and create design space
    widths = (5, 10, 15, 20, 25)
    spacing = (1, 2, 3, 4, 5)
    designSpace = product(widths, spacing)

    # define non parameter variables
    tapeAngles = (0, 90)
    tapeThickness = 0.18
    undulationRatio = 0.18
    displacements = [0.1, UNSET, UNSET]

    # run simulations
    stiffnessData = {}
    for study in designSpace:
        tapeWidth = study[0]
        tapeSpacing = study[1]
        studyStiffness = interlacedRVEImplicit.main(tapeAngles, tapeWidth,
                                                    tapeSpacing, tapeThickness,
                                                    undulationRatio, 
                                                    displacements)
        stiffnessData[study] = studyStiffness

    reportFilePath = 'C:\\Workspace\\3D_RVE\\Parametric\\PS_ResultsSummary.csv'
    with open(reportFilePath, 'wb') as csvfile:
        headings = ['tapeWidth', 'tapeSpacing', 'E11']
        csvWriter = csv.writer(csvfile, delimiter=',',
                               quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvWriter.writerow(headings)
        for key, value in stiffnessData.iteritems():
            w = key[0]
            s = key[1]
            csvWriter.writerow([w, s, value])


if __name__ == '__main__':
    parametricStudy()