from sys import path
githubPath = 'C:\\Users\\rutge\\Documents\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\mesoscale')
from abaqus import *
from abaqusConstants import *
from math import radians
import analytic_stiffness


def tapeElastic(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                G23, density):
    mat = model.Material(name='Tape-Elastic')
    mat.Elastic(type=ENGINEERING_CONSTANTS,
                table=((E11, E22, E33, nu12, nu13, nu23, G12, G13, G23), ))
    mat.Density(table=((density, ), ))
    section = model.HomogeneousSolidSection(name='Tape-Elastic',
                                            material='Tape-Elastic')


def tapeDamage(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
               G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus, G1Minus,
               G2Plus, G2Minus, G6):
    mat = model.Material(name='Tape-Damage')
    mat.Density(table=((density, ), ))
    mat.Depvar(deleteVar=20, n=24)
    mat.UserMaterial(mechanicalConstants=(E11, E22, E33, nu12, nu13, nu23,
                                          G12, G13, G23, Xt, Xc, Yt, Yc, Sl,
                                          alpha0, G1Plus, G1Minus, G2Plus,
                                          G2Minus, G6))
    section = model.HomogeneousSolidSection(name='Tape-Damage',
                                            material='Tape-Damage')


def undulationElastic(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                      G23, density, tapeAngles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    interface_angle_combos = [(0, ang) for ang in tapeAngles if ang != 0]
    n = 200  # for analytic model discretization
    for combo in interface_angle_combos:
        material_name = 'Undulation-Elastic-{}'.format(combo)
        theta_1 = radians(combo[0])  # in-plane angle of bottom ply
        theta_2 = radians(combo[1])  # in-plane angle of undulating ply
        c_lamina = analytic_stiffness.CFromConstants(
            E11, E22, nu12, nu23, G12, G23)
        v_fraction, a = analytic_stiffness.undulationGeometry(cpt, uw, n)
        c_isostrain, c_isostress = analytic_stiffness.determineStiffness(
            c_lamina, c_lamina, a, theta_1, theta_2, v_fraction, n, 1)
        eng_constants = analytic_stiffness.engineeringConstants(c_isostrain)
        mat = model.Material(name=material_name)
        mat.Elastic(type=ENGINEERING_CONSTANTS, table=(eng_constants, ))
        mat.Density(table=((density, ), ))
        section = model.HomogeneousSolidSection(name=material_name,
                                                material=material_name)


def undulationElasticResin(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                           G23, density, tapeAngles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    material_name = 'Undulation-Elastic-Resin'
    n = 200
    theta_1 = radians(0)  # in-plane angle of bottom ply
    theta_2 = radians(0)  # in-plane angle of undulating ply
    c_lamina = analytic_stiffness.CFromConstants(
        E11, E22, nu12, nu23, G12, G23)
    Clamina_resin = analytic_stiffness.CIsotropic(7.47, 0.32)
    v_fraction, a = analytic_stiffness.undulationGeometry(cpt, uw, n)
    c_isostrain, c_isostress = analytic_stiffness.determineStiffness(
        Clamina_resin, c_lamina, a, theta_1, theta_2, v_fraction, n, 2)
    eng_constants = analytic_stiffness.engineeringConstants(c_isostrain)
    material_properties = eng_constants
    mat = model.Material(name=material_name)
    mat.Elastic(type=ENGINEERING_CONSTANTS, table=(material_properties, ))
    mat.Density(table=((density, ), ))
    section = model.HomogeneousSolidSection(name=material_name,
                                            material=material_name)


def undulationDamage(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                     G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus, G1Minus,
                     G2Plus, G2Minus, G6, tapeAngles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    interface_angle_combos = [(0, ang) for ang in tapeAngles if ang != 0]
    n = 200  # for analytic model discretization
    for combo in interface_angle_combos:
        material_name = 'Undulation-Damage-{}'.format(combo)
        theta_1 = radians(combo[0])
        theta_2 = radians(combo[1])
        c_lamina = analytic_stiffness.CFromConstants(
            E11, E22, nu12, nu23, G12, G23)
        v_fraction, a = analytic_stiffness.undulationGeometry(cpt, uw, n)
        c_isostrain, c_isostress = analytic_stiffness.determineStiffness(
            c_lamina, c_lamina, a, theta_1, theta_2, v_fraction, n, 1)
        eng_constants = analytic_stiffness.engineeringConstants(c_isostrain)
        material_properties = eng_constants + [Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                                   G1Minus, G2Plus, G2Minus, G6]
        mat = model.Material(name=material_name)
        mat.Depvar(deleteVar=20, n=24)
        mat.UserMaterial(mechanicalConstants=material_properties)
        mat.Density(table=((density, ), ))
        section = model.HomogeneousSolidSection(
            name=material_name, material=material_name)


def undulationDamageResin(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                          G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                          G1Minus, G2Plus, G2Minus, G6, tapeAngles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    n = 200  # for analytic model discretization
    material_name = 'Undulation-Damage-Resin'
    theta_1 = radians(0)  # in-plane angle of bottom ply
    theta_2 = radians(0)  # in-plane angle of undulating ply
    c_lamina = analytic_stiffness.CFromConstants(
        E11, E22, nu12, nu23, G12, G23)
    Clamina_resin = analytic_stiffness.CIsotropic(7.47, 0.32)
    v_fraction, a = analytic_stiffness.undulationGeometry(cpt, uw, n)
    c_isostrain, c_isostress = analytic_stiffness.determineStiffness(
        Clamina_resin, c_lamina, a, theta_1, theta_2, v_fraction, n, 2)
    eng_constants = analytic_stiffness.engineeringConstants(c_isostrain)
    material_properties = eng_constants + [Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                               G1Minus, G2Plus, G2Minus, G6]
    mat = model.Material(name=material_name)
    mat.Depvar(deleteVar=20, n=24)
    mat.UserMaterial(mechanicalConstants=material_properties)
    mat.Density(table=((density, ), ))
    section = model.HomogeneousSolidSection(
        name=material_name, material=material_name)


def resinElastic(model):
    '''VTC401 Resin'''
    mat = model.Material(name='Resin-Elastic')
    mat.Elastic(type=ENGINEERING_CONSTANTS,
                table=((7.47, 7.47, 7.47, 0.32, 0.32, 0.32, 3.5, 3.5, 3.5), ))
    mat.Density(table=((1.1e-06, ), ))
    section = model.HomogeneousSolidSection(name='Resin-Elastic',
                                            material='Resin-Elastic')


if __name__ == '__main__':
    # Material input data
    E11 = 146.8
    E22 = E33 = 11.6
    nu12 = nu13 = 0.289
    nu23 = 0.298
    G12 = G13 = 6.47
    G23 = 4.38
    Xt = 2.61
    Xc = 1.759
    Yt = 0.055
    Yc = 0.285
    Sl = 0.105
    alpha0 = 53.0
    G1Plus = 0.1
    G1Minus = 0.1
    G2Plus = 0.00075
    G2Minus = 0.0025
    G6 = 0.0035
    density = 1.59e-06
    tapeAngles = (0, 90)
    cpt = 0.2
    uw = 1.0

    mdb.Model(name='Material Test', modelType=STANDARD_EXPLICIT)
    testModel = mdb.models['Material Test']
    tapeElastic(testModel, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                density)
    tapeDamage(testModel, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
               density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus, G1Minus, G2Plus,
               G2Minus, G6)
    undulationDamage(testModel, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                     density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus, G1Minus,
                     G2Plus, G2Minus, G6, tapeAngles, cpt, uw)
    undulationElastic(testModel, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                      density, tapeAngles, cpt, uw)
    resinElastic(testModel)