from sys import path
github_path = 'C:\\GitHub'
path.append(github_path + '\\interlaced_model_creation\\mesoscale')
path.append('C:\\Python27\\Lib\\site-packages')
from abaqus import *
from abaqusConstants import *
from math import radians, pi
from itertools import permutations
import analytic_stiffness


def tape_elastic(model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                 G23, density):
    mat = model.Material(name='Tape-Elastic')
    mat.Elastic(type=ENGINEERING_CONSTANTS,
                table=((E11, E22, E33, nu12, nu13, nu23, G12, G13, G23), ))
    mat.Density(table=((density, ), ))
    model.HomogeneousSolidSection(
        name='Tape-Elastic', material='Tape-Elastic')


def tape_damage(model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23, Xt, Xpo,
                Xc, Yt, Yc, Sl, St, alpha0, etaL, GL1Plus, GE1Plus, G1Minus,
                G2Plus, G6, density):
    mat = model.Material(name='Tape-Damage')
    mat.Density(table=((density, ), ))
    mat.Depvar(deleteVar=27, n=27)
    mat.UserMaterial(
        mechanicalConstants=(E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                             Xt, Xpo, Xc, Yt, Yc, Sl, St, alpha0, etaL,
                             GL1Plus, GE1Plus, G1Minus, G2Plus, G6))
    model.HomogeneousSolidSection(
        name='Tape-Damage', material='Tape-Damage')


def undulation_elastic(model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                       density, tape_angles, cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    resin_angles_1 = [('Resin', ang) for ang in tape_angles]
    resin_angles_2 = [(ang, 'Resin') for ang in tape_angles]
    resin_angle_combos = resin_angles_1 + resin_angles_2
    interface_angle_combos = [(0, ang) for ang in tape_angles if ang != 0]
    angle_combos = resin_angle_combos + interface_angle_combos
    n = 200  # for analytic model discretization
    # Calculate geometric parameters (constant for all angle combinations)
    v_fraction, a = analytic_stiffness.undulation_geometry(cpt, uw, n)
    # Calculate lamina stiffness matrix from engineering constants
    c_lamina = analytic_stiffness.c_from_constants(
        E11, E22, nu12, nu23, G12, G23)
    # Calculate resin region stiffness matric from engineering constants
    c_resin = analytic_stiffness.c_isotropic(7.47, 0.32)
    for combo in angle_combos:
        # Primary
        if combo[0] == 'Resin':
            material_name = "Undulation-Elastic-(Resin, {})".format(combo[1])
            c_isostrain, c_isostress = analytic_stiffness.calculate_stiffness(
                c_lamina, c_resin, a, radians(combo[1]), 0.0, v_fraction, n)
        elif combo[1] == 'Resin':
            material_name = "Undulation-Elastic-({}, Resin)".format(combo[0])
            c_isostrain, c_isostress = analytic_stiffness.calculate_stiffness(
                c_resin, c_lamina, a, 0.0, radians(combo[0]), v_fraction, n)
        else:
            material_name = "Undulation-Elastic-{}".format(combo)
            c_isostrain, c_isostress = analytic_stiffness.calculate_stiffness(
                c_lamina, c_lamina, a, radians(combo[1]), radians(combo[0]),
                v_fraction, n)
        eng_constants = analytic_stiffness.engineering_constants(c_isostrain)
        mat = model.Material(name=material_name)
        mat.Elastic(type=ENGINEERING_CONSTANTS, table=(eng_constants, ))
        mat.Density(table=((density, ), ))
        model.HomogeneousSolidSection(
            name=material_name, material=material_name)


def undulation_damage(model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                      Xt, Xpo, Xc, Yt, Yc, Sl, St, alpha0, etaL, GL1Plus,
                      GE1Plus, G1Minus, G2Plus, G6, density, tape_angles,
                      cpt, uw):
    '''
    Stiffness of undulation regions determined using analytical model
    adapted from Chou et al. 1972
    '''
    resin_angles_1 = [('Resin', ang) for ang in tape_angles]
    resin_angles_2 = [(ang, 'Resin') for ang in tape_angles]
    resin_angle_combos = resin_angles_1 + resin_angles_2
    interface_angle_combos = [(0, ang) for ang in tape_angles if ang != 0]
    angle_combos = resin_angle_combos + interface_angle_combos
    n = 200  # for analytic model discretization
    # Calculate geometric parameters (constant for all angle combinations)
    v_fraction, a = analytic_stiffness.undulation_geometry(cpt, uw, n)
    # Calculate lamina stiffness matrix from engineering constants
    c_lamina = analytic_stiffness.c_from_constants(
        E11, E22, nu12, nu23, G12, G23)
    # Calculate resin region stiffness matric from engineering constants
    c_resin = analytic_stiffness.c_isotropic(7.47, 0.32)
    for combo in angle_combos:
        # Primary
        if combo[0] == 'Resin':
            material_name = "Undulation-Damage-(Resin, {})".format(combo[1])
            c_isostrain, c_isostress = analytic_stiffness.calculate_stiffness(
                c_lamina, c_resin, a, radians(combo[1]), 0.0, v_fraction, n)
        elif combo[1] == 'Resin':
            material_name = "Undulation-Damage-({}, Resin)".format(combo[0])
            c_isostrain, c_isostress = analytic_stiffness.calculate_stiffness(
                c_resin, c_lamina, a, 0.0, radians(combo[0]), v_fraction, n)
        else:
            material_name = "Undulation-Damage-{}".format(combo)
            c_isostrain, c_isostress = analytic_stiffness.calculate_stiffness(
                c_lamina, c_lamina, a, radians(combo[1]), radians(combo[0]),
                v_fraction, n)
        eng_constants = analytic_stiffness.engineering_constants(c_isostrain)
        damage_props = [
            Xt, Xpo, Xc, Yt, Yc, Sl, St, alpha0, etaL, GL1Plus, GE1Plus,
            G1Minus, G2Plus, G6]
        material_properties = eng_constants + damage_props
        mat = model.Material(name=material_name)
        mat.Depvar(deleteVar=27, n=27)
        mat.UserMaterial(mechanicalConstants=material_properties)
        mat.Density(table=((density, ), ))
        model.HomogeneousSolidSection(
            name=material_name, material=material_name)


def resin_elastic(model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                  density):
    '''VTC401 Resin'''
    fvf = pi / 4.0
    E_m = 7.47
    nu_m = 0.32
    G_m = 3.5
    rho_m = 1.1e-06
    E_l = rule_of_mixtures(E11, E_m, fvf)
    E_t = rule_of_mixtures(E22, E_m, fvf)
    nu_l = rule_of_mixtures(nu12, nu_m, fvf)
    nu_t = rule_of_mixtures(nu23, nu_m, fvf)
    G_l = rule_of_mixtures(G12, G_m, fvf)
    G_t = rule_of_mixtures(G23, G_m, fvf)
    rho_c = rule_of_mixtures(density, rho_m, fvf)
    mat = model.Material(name='Resin-Tape-Elastic')
    eng_constants = [E_l, E_t, E_t, nu_l, nu_l, nu_t, G_l, G_l, G_t]
    mat.Elastic(type=ENGINEERING_CONSTANTS, table=(eng_constants, ))
    mat.Density(table=((rho_c, ), ))
    model.HomogeneousSolidSection(
        name='Resin-Tape-Elastic', material='Resin-Tape-Elastic')


def create_interaction_properties(model, E33, G12, GIc, GIIc, mesh_size, cpt):
    '''
    This method creates cohesive interaction properties.
    Args:
        E33 (float): transverse stiffness
        G12 (float): in plane shear modulus
        GIc (float): mode 1 fracture toughness
        GIIc (float): mode 2 fracture toughness
        mesh_size (float): size of elements in FE mesh.
        cpt (float): cured ply thickness
    Returns:
        None
    '''
    alpha = 50  # from Turon (2007)
    Ne = 5  # number of elements in the cohesive length

    # calculate cohesive zone properties according to Turon (2007)
    K1 = (alpha * E33) / cpt
    K2 = (alpha * G12) / cpt
    tau1 = ((9.0 * pi * E33 * GIc) / (32.0 * Ne * mesh_size))**0.5
    tau2 = ((9.0 * pi * E33 * GIIc) / (32.0 * Ne * mesh_size))**0.5

    model.ContactProperty('Tangential')
    model.interactionProperties['Tangential'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC,
        slipRateDependency=OFF, pressureDependency=OFF,
        temperatureDependency=OFF, dependencies=0, table=((0.15, ), ),
        shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005,
        elasticSlipStiffness=None)
    model.interactionProperties['Tangential'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON,
        constraintEnforcementMethod=DEFAULT)
    model.ContactProperty('Cohesive')
    model.interactionProperties['Cohesive'].CohesiveBehavior(
        defaultPenalties=OFF, table=((K1, K2, K2), ))
    model.interactionProperties['Cohesive'].Damage(
        criterion=QUAD_TRACTION, initTable=((tau1, tau2, tau2), ),
        useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
        mixedModeType=BK, exponent=1.75, evolTable=((GIc, GIIc, GIIc), ))


def rule_of_mixtures(prop_1, prop_2, vf):
    return vf * prop_1 + (1.0 - vf) * prop_2


if __name__ == '__main__':
    # Material input data
    E11 = 116.6
    E22 = E33 = 7.231
    nu12 = nu13 = 0.339
    nu23 = 0.374
    G12 = G13 = 3.268
    G23 = 2.632
    Xt = 2.180
    Xpo = 0.218  # estimate based on Maimi
    Xc = 0.811
    Yt = 0.131
    Yc = 0.185
    Sl = 0.122
    St = 0.070
    etaL = 0.5
    alpha0 = 53.0
    GL1Plus = 0.1034  # Tan
    GE1Plus = 0.0296  # Estimate from Tan and Maimi (same ratio)
    G1Minus = 0.095  # Furtado
    G2Plus = 0.00038
    G6 = 0.00162
    density = 1.59e-06
    tape_angles = (0, 90)
    cpt = 0.18
    uw = 1.0

    mdb.Model(name='Material Test', modelType=STANDARD_EXPLICIT)
    test_model = mdb.models['Material Test']
    tape_elastic(
        test_model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23, density)
    tape_damage(
        test_model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23, Xt, Xpo, Xc,
        Yt, Yc, Sl, St, alpha0, etaL, GL1Plus, GE1Plus, G1Minus, G2Plus, G6,
        density)
    undulation_damage(
        test_model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23, Xt, Xpo, Xc,
        Yt, Yc, Sl, St, alpha0, etaL, GL1Plus, GE1Plus, G1Minus, G2Plus, G6,
        density, tape_angles, cpt, uw)
    undulation_elastic(
        test_model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23, density,
        tape_angles, cpt, uw)
    resin_elastic(
        test_model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23, density)
    create_interaction_properties(test_model, E33, G12, G2Plus, G6, 0.5, cpt)
