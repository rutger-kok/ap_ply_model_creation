from sys import path
githubPath = 'C:\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\editing\\3D\\ELSA_PS')
import elsaSpecimenParametric as esp
from itertools import product
from abaqusConstants import *
import csv


def parametricStudy():
    # define parameters and create design space
    Lg = (30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0)
    wo = (65.0, )
    r1 = (15.0, 20.0, 25.0, 35.0)
    designSpace = product(Lg, wo, r1)

    reportFilePath = 'C:\\Workspace\\ELSA\\PS_Summary.csv'
    headings = ['Lg', 'wo', 'r1', 'r2', 'maxstress', 'sif']

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
                Lg_, wo_, r1_ = study[0], study[1], study[2]
                xdist = (D - Lg_) / 2.0
                if xdist > r1_:
                    maxstress, sif, r2 = esp.main(t, n, Lo, Lg_, D, wo_, wg, r1_)
                    csvWriter.writerow([Lg_, wo_, r1_, r2, maxstress, sif])
                else:
                    continue
            # except:
            #     csvWriter.writerow([Lg_, wo_, r1_, 'NaN', 'NaN', 'NaN'])
            #     continue


if __name__ == '__main__':
    parametricStudy()