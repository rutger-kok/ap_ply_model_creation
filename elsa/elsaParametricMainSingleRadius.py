from sys import path
githubPath = 'C:\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\editing\\3D\\ELSA_PS')
import elsaSpecimenParametricSingleRadius as esp
from itertools import product
from abaqusConstants import *
import csv


def parametricStudy():
    # define parameters and create design space
    Lg = (30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0)
    wo = (65.0, )
    designSpace = product(Lg, wo)

    reportFilePath = 'C:\\Workspace\\ELSA\\PS_Summary_SingleRadius.csv'
    headings = ['Lg', 'wo', 'r', 'maxstress', 'sif']

    # fixed parameters
    Lo = 300.0
    wg = 40.0
    D = 100.0
    t = 0.2
    n = 20

    # run simulations
    with open(reportFilePath, 'ab') as csvfile:
        csvWriter = csv.writer(csvfile, delimiter=',',
                               quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvWriter.writerow(headings)                       
        for study in designSpace:
            # try:
                Lg_, wo_ = study[0], study[1]
                maxstress, sif, r = esp.main(t, n, Lo, Lg_, D, wo_, wg)
                csvWriter.writerow([Lg_, wo_, r, maxstress, sif])
            # except:
            #     csvWriter.writerow([Lg_, wo_, 'Nan', 'NaN', 'NaN'])
            #     continue


if __name__ == '__main__':
    parametricStudy()