'''
This module is part of a library used to generate AP-PLY composite
laminate geometries in Abaqus Explicit.
Copyright (C) 2022  Rutger Kok

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
USA
'''

import sys
# change these paths to point to your local Python installation package
# libraries and the AP-PLY model creation library.
sys.path.append('C:\\Python27\\Lib\\site-packages')
sys.path.append('C:\\Github\\interlaced_model_creation')
from abaqus import *
from abaqusConstants import *
import numpy as np
from numpy import cos, radians
from itertools import permutations


def tape_damage(model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                Xt, Xc, Yt, Yc, Sl, St, G1Plus, G1Minus, G2Plus, G6,
                density):
    mat = model.Material(name='Tape-Damage')
    mat.Density(table=((density, ), ))
    mat.Depvar(deleteVar=46, n=52)
    mat.UserMaterial(
        mechanicalConstants=(E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                             Xt, Xc, Yt, Yc, Yt, Yc, Sl, Sl, St,
                             G1Plus, G1Minus, G2Plus, G6,
                             0.0, 0.0, 0.0, 0.0, 0.0))
    model.HomogeneousSolidSection(
        name='Tape-Damage', material='Tape-Damage')


def undulation_damage(model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                      Xt, Xc, Yt, Yc, Sl, St, G1Plus, G1Minus,
                      G2Plus, G6, density, tape_angles, cpt, uw):
    # define the possible angle combinations for undulations
    angle_combos = permutations(tape_angles, 2)
    for combo in angle_combos:
        material_name = "Undulation-Damage-{}".format(combo)
        mat = model.Material(name=material_name)
        mat.Density(table=((density, ), ))
        mat.Depvar(deleteVar=46, n=52)
        mat.UserMaterial(
            mechanicalConstants=(
                E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                Xt, Xc, Yt, Yc, Yt, Yc, Sl, Sl, St,
                G1Plus, G1Minus, G2Plus, G6,
                1.0, combo[0], combo[1], cpt, uw * 2.0))
        model.HomogeneousSolidSection(
            name=material_name, material=material_name)
    resin_undulation_name = "Undulation-Damage-(0, 'Resin')"
    mat = model.Material(name=resin_undulation_name)
    mat.Density(table=((density, ), ))
    mat.Depvar(deleteVar=46, n=52)
    mat.UserMaterial(
        mechanicalConstants=(
            E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            Xt, Xc, Yt, Yc, Yt, Yc, Sl, Sl, St,
            G1Plus, G1Minus, G2Plus, G6,
            1.0, 0.0, 360.0, cpt, uw * 2.0))
    model.HomogeneousSolidSection(
            name=resin_undulation_name, material=resin_undulation_name)


def resin_rich_damage(model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                      Xt, Xc, Yt, Yc, Sl, St, G1Plus, G1Minus,
                      G2Plus, G6, density, tape_angles, cpt, uw):
    # define the possible angle combinations for undulations
    material_name = "Resin-Rich-Damage"
    mat = model.Material(name=material_name)
    mat.Density(table=((density, ), ))
    mat.Depvar(deleteVar=46, n=52)
    mat.UserMaterial(
        mechanicalConstants=(
            E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            Xt, Xc, Yt, Yc, Yt, Yc, Sl, Sl, St,
            G1Plus, G1Minus, G2Plus, G6,
            2.0, 0.0, 360.0, cpt, uw * 2.0))
    model.HomogeneousSolidSection(
        name=material_name, material=material_name)


def tape_elastic(model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
                 density):
    mat = model.Material(name='Tape-Elastic')
    mat.Elastic(type=ENGINEERING_CONSTANTS,
                table=((E11, E22, E33, nu12, nu13, nu23, G12, G13, G23), ))
    mat.Density(table=((density, ), ))
    model.HomogeneousSolidSection(
        name='Tape-Elastic', material='Tape-Elastic')


def create_interaction_properties(model, tau_N, tau_S, E33, G12, GIc, GIIc,
                                  cpt):
    '''
    This method creates cohesive interaction properties.

    Args:
        tau_N (float): normal strength
        tau_S (float): shear strength
        E33 (float): transverse stiffness
        G12 (float): in plane shear modulus
        GIc (float): mode 1 fracture toughness
        GIIc (float): mode 2 fracture toughness
        cpt (float): cured ply thickness

    Returns:
        None
    '''
    alpha = 50  # from Turon (2007)

    # calculate cohesive zone properties according to Turon (2007)
    K1 = (alpha * E33) / cpt
    K2 = (alpha * G12) / cpt

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
        criterion=QUAD_TRACTION, initTable=((tau_N, tau_S, tau_S), ),
        useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
        mixedModeType=BK, exponent=1.75, evolTable=((GIc, GIIc, GIIc), ))


if __name__ == '__main__':
    # Material input data (test)
    E_11 = 124.35
    E_22 = E_33 = 7.231
    nu_12 = nu_13 = 0.339
    nu_23 = 0.374
    G_12 = G_13 = 3.268
    G_23 = 2.632
    X_t = 2.550
    X_c = 1.102
    Y_t = 0.131  # in-situ
    Y_c = 0.184
    S_l = 0.122  # in-situ
    S_t = 0.0827  # in-situ 0.0320
    G_1Plus = 0.133  # Tan
    G_1Minus = 0.095  # Furtado
    G_2Plus = 0.00038
    G_6 = 0.00162
    rho = 1.59e-06

    ply_angles = (0, 90)
    thickness = 0.205
    L_u = 1.0

    mdb.Model(name='Material Test', modelType=STANDARD_EXPLICIT)
    test_model = mdb.models['Material Test']
    tape_elastic(
        test_model, E_11, E_22, E_33, nu_12, nu_13, nu_23, G_12, G_13, G_23,
        rho)
    tape_damage(
        test_model, E_11, E_22, E_33, nu_12, nu_13, nu_23, G_12, G_13, G_23,
        X_t, X_c, Y_t, Y_c, S_l, S_t, G_1Plus, G_1Minus,
        G_2Plus, G_6, rho)
    undulation_damage(
        test_model, E_11, E_22, E_33, nu_12, nu_13, nu_23, G_12, G_13, G_23,
        X_t, X_c, Y_t, Y_c, S_l, S_t, G_1Plus, G_1Minus,
        G_2Plus, G_6, rho, ply_angles, thickness, L_u)
    create_interaction_properties(
        test_model, Y_t, S_l, E_33, G_12, G_2Plus, G_6, thickness)
