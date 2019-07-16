import numpy as np
import math

def rotateLaminaStiffness(C_lam, O):
    # define the transformation matrix for in-plane rotation
    Orad = math.radians(O)
    m = math.cos(Orad)
    n = math.sin(Orad)

    T_ip = np.array(
        [(m**2, n**2, 0.0, 0.0, 0.0, 2*m*n),
        (n**2, m**2, 0.0, 0.0, 0.0, -2*m*n),
        (0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, m, -n, 0.0),
        (0.0, 0.0, 0.0, n, m, 0.0),
        (-m*n, m*n, 0.0, 0.0, 0.0, (m**2)-(n**2))])
    T_ip_inv = np.linalg.inv(T_ip) # invert in-plane rotation matrix

    # create the Reuter matrix
    R = np.identity(6) 
    for g in range(3,6):  # zero based indexing!
        R[g, g] = 2.0
    R_inv = np.linalg.inv(R)

    ip_rotation = np.dot(np.dot(np.dot(np.dot(
            T_ip_inv, C_lam), R), T_ip), R_inv)
    return ip_rotation

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

def rotateMatProps(E_11, E_22, nu_12, nu_23, G_12, G_23, angle):
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

    S_lamina = np.array(
            [(S11, S12, S13, 0.0, 0.0, 0.0),
            (S21, S22, S23, 0.0, 0.0, 0.0),
            (S31, S32, S33, 0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, S44, 0.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, S55, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, S66)])

    C_lamina = np.linalg.inv(S_lamina)
    print C_lamina
    S_rotated = np.linalg.inv(rotateLaminaStiffness(C_lamina, angle))
    matprops_rotated = engineeringConstants(S_rotated)

    return matprops_rotated

if __name__ == '__main__':
    print 'Testing'
    print rotateMatProps(161.0, 11.4, 0.32, 0.436, 5.29, 3.98, 0.0)
    # print rotateMatProps(125.0, 8.56, 0.32, 0.23, 6.47, 4.5, 90.0)
    # print rotateMatProps(140.40, 11.61, 0.289, 0.299, 6.47, 4.38, 45)