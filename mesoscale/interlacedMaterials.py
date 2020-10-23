from sys import path
githubPath = 'C:\\Users\\rutge\\Documents\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\editing')
from abaqus import *
from abaqusConstants import *
from math import radians
import analyticStiffness as analyticStiffness


def tapeElastic(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                G23, density):
    mat = model.Material(name='Tape-Elastic')
    mat.Elastic(type=ENGINEERING_CONSTANTS,
                table=((E11, E22, E33, nu12, nu13, nu23, G12, G13, G23), ))
    mat.Density(table=((density, ), ))
    section = model.HomogeneousSolidSection(name='Tape-Elastic-Section',
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
    section = model.HomogeneousSolidSection(name='Tape-Damage-Section',
                                            material='Tape-Damage')


def undulationElastic(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                      G23, density, tapeAngles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    interfaceAngleCombos = [(0, ang) for ang in tapeAngles if ang != 0]
    n = 200  # for analytic model discretization
    for combo in interfaceAngleCombos:
        matName = 'Undulation-Elastic-{}'.format(combo)
        O1 = radians(combo[0])  # in-plane angle of bottom ply
        O2 = radians(combo[1])  # in-plane angle of undulating ply
        Clamina = analyticStiffness.CFromConstants(
            E11, E22, nu12, nu23, G12, G23)
        Vfrac, a = analyticStiffness.undulationGeometry(cpt, uw, n)
        CIsostrain, CIsostress = analyticStiffness.determineStiffness(
            Clamina, Clamina, a, O1, O2, Vfrac, n, 1)
        engConstants = analyticStiffness.engineeringConstants(CIsostrain)
        matProps = engConstants
        mat = model.Material(name=matName)
        mat.Elastic(type=ENGINEERING_CONSTANTS, table=(matProps, ))
        mat.Density(table=((density, ), ))
        section = model.HomogeneousSolidSection(name=matName + '-Section',
                                                material=matName)


def undulationElasticResin(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                           G23, density, tapeAngles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    matName = 'Undulation-Elastic-Resin'
    n = 200
    O1 = radians(0)  # in-plane angle of bottom ply
    O2 = radians(0)  # in-plane angle of undulating ply
    Clamina = analyticStiffness.CFromConstants(
        E11, E22, nu12, nu23, G12, G23)
    Clamina_resin = analyticStiffness.CIsotropic(7.47, 0.32)
    Vfrac, a = analyticStiffness.undulationGeometry(cpt, uw, n)
    CIsostrain, CIsostress = analyticStiffness.determineStiffness(
        Clamina_resin, Clamina, a, O1, O2, Vfrac, n, 2)
    engConstants = analyticStiffness.engineeringConstants(CIsostrain)
    matProps = engConstants
    mat = model.Material(name=matName)
    mat.Elastic(type=ENGINEERING_CONSTANTS, table=(matProps, ))
    mat.Density(table=((density, ), ))
    section = model.HomogeneousSolidSection(name=matName + '-Section',
                                            material=matName)


def undulationDamage(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                     G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus, G1Minus,
                     G2Plus, G2Minus, G6, tapeAngles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    interfaceAngleCombos = [(0, ang) for ang in tapeAngles if ang != 0]
    n = 200  # for analytic model discretization
    for combo in interfaceAngleCombos:
        matName = 'Undulation-Damage-{}'.format(combo)
        O1 = radians(combo[0])
        O2 = radians(combo[1])
        Clamina = analyticStiffness.CFromConstants(
            E11, E22, nu12, nu23, G12, G23)
        Vfrac, a = analyticStiffness.undulationGeometry(cpt, uw, n)
        CIsostrain, CIsostress = analyticStiffness.determineStiffness(
            Clamina, Clamina, a, O1, O2, Vfrac, n, 1)
        engConstants = analyticStiffness.engineeringConstants(CIsostrain)
        matProps = engConstants + [Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                                   G1Minus, G2Plus, G2Minus, G6]
        mat = model.Material(name=matName)
        mat.Depvar(deleteVar=20, n=24)
        mat.UserMaterial(mechanicalConstants=matProps)
        mat.Density(table=((density, ), ))
        section = model.HomogeneousSolidSection(
            name=matName + '-Section', material=matName)


def undulationDamageResin(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                          G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                          G1Minus, G2Plus, G2Minus, G6, tapeAngles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    n = 200  # for analytic model discretization
    matName = 'Undulation-Damage-Resin'
    O1 = radians(0)  # in-plane angle of bottom ply
    O2 = radians(0)  # in-plane angle of undulating ply
    Clamina = analyticStiffness.CFromConstants(
        E11, E22, nu12, nu23, G12, G23)
    Clamina_resin = analyticStiffness.CIsotropic(7.47, 0.32)
    Vfrac, a = analyticStiffness.undulationGeometry(cpt, uw, n)
    CIsostrain, CIsostress = analyticStiffness.determineStiffness(
        Clamina_resin, Clamina, a, O1, O2, Vfrac, n, 2)
    engConstants = analyticStiffness.engineeringConstants(CIsostrain)
    matProps = engConstants + [Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                               G1Minus, G2Plus, G2Minus, G6]
    mat = model.Material(name=matName)
    mat.Depvar(deleteVar=20, n=24)
    mat.UserMaterial(mechanicalConstants=matProps)
    mat.Density(table=((density, ), ))
    section = model.HomogeneousSolidSection(
        name=matName + '-Section', material=matName)


def resinElastic(model):
    '''VTC401 Resin'''
    mat = model.Material(name='Resin-Elastic')
    mat.Elastic(type=ENGINEERING_CONSTANTS,
                table=((7.47, 7.47, 7.47, 0.32, 0.32, 0.32, 3.5, 3.5, 3.5), ))
    mat.Density(table=((1.1e-06, ), ))
    section = model.HomogeneousSolidSection(name='Resin-Elastic-Section',
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