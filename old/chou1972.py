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
a tape undulation.
Inputs to the model are as follows:
Geometric:
ht = 0.8  - total height of the laminate (mm)
hf = 0.2  - height of undulating plies (mm)
hu = 0.2  - height of the undulation (mm)
Lu = 1.0  - total length of the undulation (mm)
O1 = math.radians(90)  # in-plane angle of non-undulating plies
O2 = math.radians(0)  # in-plane angle of undulating plies
Material:
Clamina = lamina stiffness matrix (can be derived from engineering constants)
Numerical:
n = 200  # number of increments
'''

def rotateLaminaStiffness(Ck, theta, alpha): # returns Ck_rot, rotated lamina C
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
    T_ip_inv = np.linalg.inv(T_ip) # invert in-plane rotation matrix

    # define the transformation matrix for out of plane rotation
    k = math.cos(alpha)
    l = math.sin(alpha)
    T_oop = np.array(
        [(k**2, 0.0, l**2, 0.0, -2*k*l, 0.0),
        (0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
        (l**2, 0.0, k**2, 0.0, 2*k*l, 0.0),
        (0.0, 0.0, 0.0, k, 0.0, l),
        (k*l, 0.0, -k*l, 0.0, (k**2)-(l**2), 0.0),
        (0.0, 0.0, 0.0, -l, 0.0, k)])
    T_oop_inv = np.linalg.inv(T_oop)

    # create the Reuter matrix
    R = np.identity(6) 
    for g in range(3,6):  # zero based indexing!
        R[g, g] = 2.0
    R_inv = np.linalg.inv(R)

    # create the rotated stiffness matrix for the lamina, accounting for both
    # in-plane and out-of-plane rotation
    ip_rotation = np.dot(np.dot(np.dot(np.dot(
            T_ip_inv, Ck), R), T_ip), R_inv)
    Ck_rot = np.dot(np.dot(np.dot(np.dot(
            T_oop_inv, ip_rotation), R), T_oop), R_inv)
    # transpose_ip = np.dot(np.dot(T_ip,C_lam),T_ip_inv)
    # transpose_oop = np.dot(np.dot(T_oop,transpose_ip),T_oop_inv)

    return Ck_rot

def engineeringConstants(Smatrix): # returns eng_const, engineering constants
    Ex = 1/Smatrix[0,0]
    Ey = 1/Smatrix[1,1]
    Ez = 1/Smatrix[2,2]
    nu_xy = -Smatrix[1,0]/Smatrix[0,0]
    nu_xz = -Smatrix[2,0]/Smatrix[0,0]
    nu_yz = -Smatrix[1,2]/Smatrix[1,1]
    Gxy = 1/Smatrix[5,5]
    Gxz = 1/Smatrix[4,4]
    Gyz = 1/Smatrix[3,3]
    eng_const = (Ex, Ey, Ez, nu_xy, nu_xz, nu_yz, Gxy, Gxz, Gyz)
    return eng_const

def undulationGeometry(h_t, h_f, h_u, L_u, n): # returns V, alpha, ratio
    # adapted for application to Abaqus; non-symmetric undulation
    x = np.linspace(0, L_u, n)  # discretize x
    ratio = h_u/L_u

    alpha = np.arctan(((math.pi*h_u)/(2.0*L_u))*np.cos((math.pi/L_u)*(x - L_u/2.0)))

    h_1 = [0.0 for qq in range(n)]
    h_2 = h_t/4.0 - h_f/2.0 + (h_u/2.0)*np.sin((math.pi/L_u)*(x - L_u/2.0))
    h_3 = h_t/4.0 + h_f/2.0 + (h_u/2.0)*np.sin((math.pi/L_u)*(x - L_u/2.0))
    h_4 = [h_t/2.0 for ww in range(n)]

    # ply volume fractions
    Vk1 = (h_4-h_3)/h_t
    Vk2 = (h_3-h_2)/h_t
    Vk3 = (h_2-h_1)/h_t
    V = [[Vk1[s], Vk2[s], Vk3[s]] for s in range(n)]

    return V, alpha, ratio

def determineStiffness(C_lamina, alpha, O_1, O_2, V, n):
    C = []
    for m in range(n):
        Ck1 = rotateLaminaStiffness(C_lamina, O_1, 0.0)
        Ck2 = rotateLaminaStiffness(C_lamina, O_2, alpha[m])
        Ck = [Ck1, Ck2, Ck1]
        Cx = np.zeros((6,6), dtype=float)
        for i in (0,1,2,5):
            for j in (0,1,2,5):
                Cx[i,j] = sum([(V[m][k]*(Ck[k][i,j]-((Ck[k][i,2]*Ck[k][2,j])/Ck[k][2,2]) + 
                    (Ck[k][i,2]/Ck[k][2,2])*(sum([((V[m][l]*Ck[l][2,j])/Ck[l][2,2]) for l in range(len(Ck))])/
                    sum([(V[m][f]/Ck[f][2,2]) for f in range(len(Ck))])))) for k in range(len(Ck))])
        for g in (3,4):
                Cx[i,g] = Cx[g,i] = 0.0
        for u in (3,4):
            for v in (3,4):
                num = 0
                den = 0
                for z in range(len(Ck)):
                    num  += (V[m][z]/np.linalg.det(np.array([[Ck[z][3,3],Ck[z][3,4]],[Ck[z][4,3],Ck[z][4,4]]])))*Ck[z][u,v]
                    for t in range(len(Ck)):
                        den += ((
                            V[m][z]/np.linalg.det(np.array([[Ck[z][3,3],Ck[z][3,4]],[Ck[z][4,3],Ck[z][4,4]]])))
                            *(V[m][t]/np.linalg.det(np.array([[Ck[t][3,3],Ck[t][3,4]],[Ck[t][4,3],Ck[t][4,4]]])))
                            *(Ck[z][3,3]*Ck[t][4,4]-Ck[z][3,4]*Ck[t][4,3]))
                Cx[u,v] = num/den
        C.append(Cx)

    # Iso-strain homogenization (averaging stiffness)
    C_isostrain = np.zeros((6,6), dtype=float)
    for ii in range(6):
        for jj in range(6):
            Ciijj = [C[d][ii,jj] for d in range(n)]
            C_isostrain[ii,jj] = sum(Ciijj)/len(Ciijj)
    # print 'Iso-strain'
    # print6x6(C_isostrain)

    # Iso-stress RVE homogenization (averaging compliance)
    S = []
    for ss in range(n):
        S.append(np.linalg.inv(C[ss]))
    S_homo = np.zeros((6,6), dtype=float)
    for oo in range(6):
        for pp in range(6):
            Soopp = [S[dd][oo,pp] for dd in range(n)]
            S_homo[oo,pp] = sum(Soopp)/len(Soopp)
    C_isostress = np.linalg.inv(S_homo)
    # print 'Iso-stress'
    # print6x6(C_isostress)
    # stiffness_isostrain = engineeringConstants(np.linalg.inv(C_isostrain))
    # stiffness_isostress = engineeringConstants(S_homo)
    
    return C_isostrain, C_isostress

def CFromConstants(E_11,E_22,nu_12,nu_23,G_12,G_23):
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
            rc1 = round(c1,2)
            rc2 = round(c2,2)
            rc3 = round(c3,2)
            rc4 = round(c4,2)
            rc5 = round(c5,2)
            rc6 = round(c6,2)
            print "{:6s}|{:6s}|{:6s}|{:6s}|{:6s}|{:6s}".format(
                str(rc1), str(rc2), str(rc3), str(rc4), str(rc5), str(rc6))
        print '_'*41

    Clamina = CFromConstants(140.4,11.6,0.289,0.298,6.47,4.38)
    n = 200  # number of increments
    ht = 0.4  # total height of the laminate
    hf = 0.2  # height of undulating plies (ply thickness * num of plies)
    hu = 0.2  # height of the undulation
    Lu = 1.0  # total length of the undulation
    bO1 = math.radians(90)  # in-plane angle of non-undulating plies
    bO2 = math.radians(0)  # in-plane angle of undulating plies

    bVfrac, ba, brat = undulationGeometry(ht, hf, hu, Lu, n)
    btemp1, btemp2 = determineStiffness(Clamina, ba, bO1, bO2, bVfrac, n)
