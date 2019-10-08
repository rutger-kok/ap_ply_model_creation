# 3D stiffness rotation
import numpy as np
import math

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

def CFromConstants(E_11,E_22,E_33,nu_12,nu_13,nu_23,G_12,G_13,G_23):
    # NOTE: this assumes transverse isotropy! (12 and 13 are the same)
    # Material Parameters
    nu_21 = nu_12*(E_22/E_11)
    nu_31 = nu_13*(E_33/E_22)
    nu_32 = nu_23*(E_33/E_22)

    # Compliance
    S11 = 1/E_11
    S22 =  1/E_22
    S33 = 1/E_33
    S12 =  -nu_12/E_11
    S13 = -nu_13/E_11
    S23 = -nu_32/E_33
    S21 = S12
    S31 = S13
    S32 = S23
    S66 = 1/G_12
    S44 = 1/G_23
    S55 = 1/G_13

    Slamina = np.array(
            [(S11, S12, S13, 0.0, 0.0, 0.0),
            (S21, S22, S23, 0.0, 0.0, 0.0),
            (S31, S32, S33, 0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, S44, 0.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, S55, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, S66)])

    Clamina = np.linalg.inv(Slamina)
    return Clamina
