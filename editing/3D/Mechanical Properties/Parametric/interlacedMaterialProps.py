from abaqus import *
from abaqusConstants import *
import analyticStiffness
from math import radians

# Material input data (VTC401)
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


def VTC401_UserMaterial(modelName, tapeAngles, tapeThickness, undulationWidth):
    activeModel = mdb.models[modelName]

    # Create Abaqus materials
    tapeMaterial = activeModel.Material(name='Tape')
    tapeMaterial.Density(table=((density, ), ))
    tapeMaterial.Depvar(deleteVar=20, n=24)
    tapeMaterial.UserMaterial(mechanicalConstants=(E11, E22, E33, nu12, nu13,
                                                   nu23, G12, G13, G23, Xt, Xc,
                                                   Yt, Yc, Sl, alpha0, G1Plus,
                                                   G1Minus, G2Plus, G2Minus,
                                                   G6))
    tapeSection = activeModel.HomogeneousSolidSection(name='Tape Section',
                                                      material='Tape')

    resinMaterial = activeModel.Material(name='Resin')
    resinMaterial.Elastic(type=ENGINEERING_CONSTANTS,
                          table=((7.47, 7.47, 7.47, 0.32, 0.32, 0.32, 3.5, 3.5,
                                  3.5), ))
    resinMaterial.Density(table=((1.1e-06, ), ))
    resinSection = activeModel.HomogeneousSolidSection(name='Resin Section',
                                                       material='Resin')

    # Determine stiffness of undulation regions using analytical model
    # Adapted from Chou et al. 1972

    interfaceAngleCombos = [(0, ang) for ang in tapeAngles if ang != 0]
    for combo in interfaceAngleCombos:
        matName = 'Undulation {}'.format(combo)
        O1 = radians(combo[1])  # in-plane angle of bottom ply
        O2 = radians(combo[0])  # in-plane angle of undulating ply
        Clamina = analyticStiffness.CFromConstants(E11, E22, nu12, nu23, G12,
                                                   G23)
        Vfrac, a = analyticStiffness.undulationGeometry(tapeThickness,
                                                        undulationWidth, n)
        CIsostrain, CIsostress = analyticStiffness.determineStiffness(
            Clamina, a, O1, O2, Vfrac, n)
        engConstants = analyticStiffness.engineeringConstants(CIsostrain)
        matProps = engConstants + [2.61, 1.759, 0.055, 0.285, 0.105, 53.0, 0.1,
                                   0.1, 0.00075, 0.0025, 0.0035]
        undMaterial = activeModel.Material(name=matName)
        undMaterial.Depvar(deleteVar=20, n=24)
        undMaterial.UserMaterial(mechanicalConstants=matProps)
        undMaterial.Density(table=((density, ), ))
        undulationSection = activeModel.HomogeneousSolidSection(
                name=matName + ' Section', material=matName)


def VTC401_Elastic(modelName, tapeAngles, tapeThickness, undulationWidth):
    activeModel = mdb.models[modelName]

    # Create Abaqus materials

    tapeMaterial = activeModel.Material(name='Tape')
    tapeMaterial.Density(table=((density, ), ))
    tapeMaterial.Elastic(type=ENGINEERING_CONSTANTS,
                         table=((E11, E22, E33, nu12, nu13, nu23, G12, G13,
                                G23), ))
    tapeSection = activeModel.HomogeneousSolidSection(name='Tape Section',
                                                      material='Tape')

    resinMaterial = activeModel.Material(name='Resin')
    resinMaterial.Elastic(type=ENGINEERING_CONSTANTS,
                          table=((7.47, 7.47, 7.47, 0.32, 0.32, 0.32, 3.5, 3.5,
                                  3.5), ))
    resinMaterial.Density(table=((1.1e-06, ), ))
    resinSection = activeModel.HomogeneousSolidSection(name='Resin Section',
                                                       material='Resin')

    # Determine stiffness of undulation regions using analytical model
    # Adapted from Chou et al. 1972

    interfaceAngleCombos = [(0, ang) for ang in tapeAngles if ang != 0]
    n = 200  # for analytic model discretization
    for combo in interfaceAngleCombos:
        matName = 'Undulation {}'.format(combo)
        O1 = radians(combo[1])  # in-plane angle of bottom ply
        O2 = radians(combo[0])  # in-plane angle of undulating ply
        Clamina = analyticStiffness.CFromConstants(E11, E22, nu12, nu23, G12,
                                                   G23)
        Vfrac, a = analyticStiffness.undulationGeometry(tapeThickness,
                                                        undulationWidth, n)
        CIsostrain, CIsostress = analyticStiffness.determineStiffness(
            Clamina, a, O1, O2, Vfrac, n)
        engConstants = analyticStiffness.engineeringConstants(CIsostrain)
        matProps = engConstants
        undMaterial = activeModel.Material(name=matName)
        undMaterial.Elastic(type=ENGINEERING_CONSTANTS, table=(matProps, ))
        undMaterial.Density(table=((density, ), ))
        undulationSection = activeModel.HomogeneousSolidSection(
                name=matName + ' Section', material=matName)