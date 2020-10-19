# ****************************************************************************
# 3D Homogenized Undulation Model adapted from Chou 1972
# Author: Rutger Kok
# Last modified: 12/08/2019
# ****************************************************************************

import math
import numpy as np

'''
This script uses the Chou(1972)/Bogetti et al. (1994) model to determine the
3D iso-stress and iso-strain homgenized stiffness of a unit cell containing
part of a tape undulation.

Inputs to the model are as follows:
Geometric:
cpt = total height of the laminate (mm)
Lu = total length of the undulation (mm)
O1 = in-plane angle of non-undulating plies
O2 = in-plane angle of undulating plies
Material:
Clamina = lamina stiffness matrix (can be derived from engineering constants)
Numerical:
n = 200  # number of increments
'''


def rotateLaminaStiffness(Ck, theta, alpha):
    # returns Ck_rot, rotated lamina C
    # define the transformation matrix for in-plane rotation
    m = math.cos(theta)
    n = math.sin(theta)
    T_ip = np.array(
                    [(m**2, n**2, 0.0, 0.0, 0.0, 2*m*n),
                     (n**2, m**2, 0.0, 0.0, 0.0, -2*m*n),
                     (0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
                     (0.0, 0.0, 0.0, m, -n, 0.0),
                     (0.0, 0.0, 0.0, n, m, 0.0),
                     (-m*n, m*n, 0.0, 0.0, 0.0, (m**2)-(n**2))])
    T_ip_inv = np.linalg.inv(T_ip)  # invert in-plane rotation matrix

    # define the transformation matrix for out-of-plane rotation
    g = math.cos(alpha)
    h = math.sin(alpha)
    T_oop = np.array(
                     [(g**2, 0.0, h**2, 0.0, -2*g*h, 0.0),
                      (0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
                      (h**2, 0.0, g**2, 0.0, 2*g*h, 0.0),
                      (0.0, 0.0, 0.0, g, 0.0, h),
                      (g*h, 0.0, -g*h, 0.0, (g**2)-(h**2), 0.0),
                      (0.0, 0.0, 0.0, -h, 0.0, g)])
    T_oop_inv = np.linalg.inv(T_oop)

    # create the Reuter matrix
    R = np.identity(6)
    for a in range(3, 6):  # zero based indexing!
        R[a, a] = 2.0
    R_inv = np.linalg.inv(R)

    # create the rotated stiffness matrix for the lamina, accounting for both
    # in-plane and out-of-plane rotation
    ip_rotation = np.dot(np.dot(np.dot(np.dot(
            T_ip_inv, Ck), R), T_ip), R_inv)
    Ck_rot = np.dot(np.dot(np.dot(np.dot(
            T_oop_inv, ip_rotation), R), T_oop), R_inv)

    return Ck_rot


def engineeringConstants(Cmatrix):  # returns eng_const, engineering constants
    Smatrix = np.linalg.inv(Cmatrix)
    Ex = 1/Smatrix[0, 0]
    Ey = 1/Smatrix[1, 1]
    Ez = 1/Smatrix[2, 2]
    nu_xy = -Smatrix[1, 0]/Smatrix[0, 0]
    nu_xz = -Smatrix[2, 0]/Smatrix[0, 0]
    nu_yz = -Smatrix[1, 2]/Smatrix[1, 1]
    Gxy = 1/Smatrix[5, 5]
    Gxz = 1/Smatrix[4, 4]
    Gyz = 1/Smatrix[3, 3]
    eng_const = [Ex, Ey, Ez, nu_xy, nu_xz, nu_yz, Gxy, Gxz, Gyz]
    return eng_const


def undulationGeometry(cpt, L_u, n):  # returns V, alpha, ratio
    # adapted for application to Abaqus; non-symmetric undulation
    x = np.linspace(1 * 10**-12, L_u, n)  # discretize x

    h_1 = [0.0 for _ in range(n)]
    h_2 = 0.5*cpt + 0.5*cpt*np.sin((math.pi/2*L_u)*(x - L_u))
    h_3 = [0.5*cpt for _ in range(n)]

    alpha = np.arctan(h_2 / x)
    alpha[np.isnan(alpha)] = 0.0  # replace NaNs with 0.0

    # ply volume fractions
    Vk1 = (h_2-h_1)/cpt
    Vk2 = (h_3-h_2)/cpt
    V = [[Vk1[s], Vk2[s]] for s in range(n)]

    return V, alpha


def determineStiffness(C_lamina1, C_lamina2, alpha, O_1, O_2, V, n):
    C = []
    for a in range(n):
        Ck1 = rotateLaminaStiffness(C_lamina1, O_1, 0.0)
        Ck2 = rotateLaminaStiffness(C_lamina2, O_2, alpha[a])
        Ck = [Ck1, Ck2]
        Cx = np.zeros((6, 6), dtype=float)
        for i in (0, 1, 2, 5):
            for j in (0, 1, 2, 5):
                Cx[i, j] = sum([(V[a][k]*(Ck[k][i, j] -
                                 ((Ck[k][i, 2]*Ck[k][2, j])/Ck[k][2, 2]) +
                                 (Ck[k][i, 2]/Ck[k][2, 2]) *
                                 (sum([((V[a][l]*Ck[l][2, j])/Ck[l][2, 2])
                                  for l in range(len(Ck))]) /
                                  sum([(V[a][f]/Ck[f][2, 2])
                                      for f in range(len(Ck))]))))
                                for k in range(len(Ck))])
        for b in (3, 4):
                Cx[i, b] = Cx[b, i] = 0.0
        for u in (3, 4):
            for v in (3, 4):
                num = 0
                den = 0
                for z in range(len(Ck)):
                    num += (V[a][z]/np.linalg.det(
                            np.array([[Ck[z][3, 3], Ck[z][3, 4]],
                                     [Ck[z][4, 3], Ck[z][4, 4]]])))*Ck[z][u, v]
                    for t in range(len(Ck)):
                        den += ((
                                 V[a][z]/np.linalg.det(
                                    np.array([[Ck[z][3, 3], Ck[z][3, 4]],
                                             [Ck[z][4, 3], Ck[z][4, 4]]])))
                                * (V[a][t]/np.linalg.det(
                                    np.array([[Ck[t][3, 3], Ck[t][3, 4]],
                                             [Ck[t][4, 3], Ck[t][4, 4]]])))
                                * (Ck[z][3, 3]*Ck[t][4, 4] -
                                   Ck[z][3, 4]*Ck[t][4, 3]))
                Cx[u, v] = num/den
        C.append(Cx)

    # Iso-strain homogenization (averaging stiffness)
    C_isostrain = np.zeros((6, 6), dtype=float)
    for x in range(6):
        for y in range(6):
            Cxy = [C[d][x, y] for d in range(n)]
            C_isostrain[x, y] = sum(Cxy)/len(Cxy)

    # Iso-stress RVE homogenization (averaging compliance)
    S = []
    for f in range(n):
        S.append(np.linalg.inv(C[f]))
    S_homo = np.zeros((6, 6), dtype=float)
    for e in range(6):
        for f in range(6):
            Sef = [S[dd][e, f] for dd in range(n)]
            S_homo[e, f] = sum(Sef)/len(Sef)
    C_isostress = np.linalg.inv(S_homo)

    return C_isostrain, C_isostress


def CFromConstants(E_11, E_22, nu_12, nu_23, G_12, G_23):
    # NOTE: this assumes transverse isotropy! (12 and 13 are the same)
    # Material Parameters
    # VTC401 T700 props as estimated using Autodesk Composite
    E_33 = E_22
    nu_13 = nu_12
    nu_21 = nu_12*(E_22/E_11)
    nu_31 = nu_13*(E_33/E_22)
    nu_32 = nu_23*(E_33/E_22)

    # Compliance
    S11 = 1/E_11
    S22 = S33 = 1/E_22
    S12 = S13 = -nu_12/E_11
    S23 = -nu_32/E_33
    S21 = S12
    S31 = S13
    S32 = S23
    S66 = S55 = 1/G_12
    S44 = 1/G_23

    Slamina = np.array(
                       [(S11, S12, S13, 0.0, 0.0, 0.0),
                        (S21, S22, S23, 0.0, 0.0, 0.0),
                        (S31, S32, S33, 0.0, 0.0, 0.0),
                        (0.0, 0.0, 0.0, S44, 0.0, 0.0),
                        (0.0, 0.0, 0.0, 0.0, S55, 0.0),
                        (0.0, 0.0, 0.0, 0.0, 0.0, S66)])

    Clamina = np.linalg.inv(Slamina)
    return Clamina


# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    def print6x6(matrix):  # function to print a more clear 6x6 matrix
        print '_'*41
        for c1, c2, c3, c4, c5, c6 in matrix:
            rc1 = round(c1, 2)
            rc2 = round(c2, 2)
            rc3 = round(c3, 2)
            rc4 = round(c4, 2)
            rc5 = round(c5, 2)
            rc6 = round(c6, 2)
            print "{:6s}|{:6s}|{:6s}|{:6s}|{:6s}|{:6s}".format(
                str(rc1), str(rc2), str(rc3), str(rc4), str(rc5), str(rc6))
        print '_'*41

    Clamina = CFromConstants(140.4, 11.6, 0.289, 0.298, 6.47, 4.38)
    n = 200  # number of increments
    cpt = 0.2
    L_u = 1.0  # total length of the undulation
    O1 = math.radians(90)  # in-plane angle of non-undulating plies
    O2 = math.radians(0)  # in-plane angle of undulating plies

    Vfrac, a = undulationGeometry(cpt, L_u, n)
    temp1, temp2 = determineStiffness(Clamina, Clamina, a, O1, O2, Vfrac, n)
    print6x6(temp1)