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
sys.path.append('C:\\Github\\ap_ply_model_creation')
from abaqus import *
from abaqusConstants import *
from shapely.geometry import Polygon
import ap_ply_materials as mats
import tape_placement as tp
from abaqus_model import Abaqus_Model


class AP_PLY_Model(Abaqus_Model):
    def __init__(self, model_name, dir_name=None):
        Abaqus_Model.__init__(self, model_name, dir_name)

    def set_specimen_parameters(self, t_angles, t_widths, t_spacing,
                                t_thickness, u_ratio, l_plies, x_shift=0.0,
                                y_shift=0.0):
        '''
        Set AP-PLY specimen parameters. Note t_xxx indicates a tape
        property, u_xxx indicates an undulation property, and l_xxx indicates
        a laminate property

        Args:
            t_angles (list): angles of tapes in the AP-PLY laminate
            t_widths (float): width of the tapes in the laminate
            t_spacing (int): number of tape width gaps between tapes placed
                in a single pass
            t_thickness (float): thickness of the tapes in the laminate
            u_ratio (float): ratio of the height of an undulation to its length
            l_plies (int): number of plies in the laminate.

        Returns:
            None
        '''
        self.t_angles = t_angles
        self.t_widths = [t_widths] * len(self.t_angles)
        self.t_spacing = t_spacing
        self.t_thickness = t_thickness
        self.u_ratio = u_ratio
        self.l_plies = l_plies  # number of plies in laminate
        self.sequences = int((l_plies / int(len(self.t_angles))) / 2)
        self.u_width = (self.t_thickness / self.u_ratio) / 2.0
        self.l_thickness = l_plies * self.t_thickness
        self.x_shift = x_shift
        self.y_shift = y_shift

    def create_materials(self, material_name=None):
        '''
        Create materials required for AP-PLY models. Makes use of the
        ap_ply_materials module.

        Args:
            None

        Returns:
            None
        '''

        if material_name is None:  # default to SHD Composites VTC401
            E11 = 124.35
            E22 = E33 = 7.231
            nu12 = nu13 = 0.339
            nu23 = 0.374
            G12 = G13 = 3.268
            G23 = 2.632
            Xt = 2.550
            Xc = 1.102
            # Yt = 0.04044
            Yt = 0.131
            # Yc = 0.0841
            # Sl = 0.0554
            Sl = 0.122
            # Yt = 0.131  # in-situ
            Yc = 0.184
            # Sl = 0.122  # in-situ
            St = 0.08273
            G1Plus = 0.133  # Tan
            G1Minus = 0.095  # Furtado
            G2Plus = 0.00038
            G6 = 0.00162
            density = 1.59e-06
        elif material_name == 'HiTapeUD210':  # Hexcel HiTape UD210
            E11 = 156.0
            E22 = E33 = 7.24
            nu12 = nu13 = 0.309
            nu23 = 0.401
            G12 = G13 = 4.65
            G23 = 3.10
            Xt = 3.135
            Xc = 1.760
            Yt = 0.076
            Sl = 0.101
            Yc = 0.180  # chamis from Herraez braid paper
            St = 0.101
            G1Plus = 0.133  # Tan
            G1Minus = 0.095  # Furtado
            G2Plus = 0.00078
            G6 = 0.00273  # 1.2 x Israr values (T700/M21)
            density = 1.59e-06
        else:
            print 'Not a valid material name'

        mats.tape_elastic(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            density)
        mats.tape_damage(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            Xt, Xc, Yt, Yc, Sl, St, G1Plus, G1Minus, G2Plus, G6,
            density)
        mats.undulation_damage(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            Xt, Xc, Yt, Yc, Sl, St, G1Plus, G1Minus, G2Plus, G6,
            density, self.t_angles, self.t_thickness, self.u_width)
        mats.resin_rich_damage(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            Xt, Xc, Yt, Yc, Sl, St, G1Plus, G1Minus, G2Plus, G6,
            density, self.t_angles, self.t_thickness, self.u_width)
        mats.create_interaction_properties(
            self.model, Yt, Sl, E33, G12, G2Plus, G6, self.t_thickness)

    def create_tape_paths(self, specimen_size):
        '''
        Uses the shapely package to define the AP-PLY specimen geometry.

        Args:
            specimen_size (Shapely Polygon): Polygon defining the dimensions of
                the specimen to be modeled.

        Returns:
            tape_paths (list): List of tape/tow connectivity for each tape.
            trimmed_grid (dict): Nested dictionary defining the grid of
                tape/undulation objects.
        '''
        self.specimen_size = specimen_size
        x, y = zip(*self.specimen_size.exterior.coords)
        x_max = max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)
        x_mid = (x_max + x_min) / 2.0
        y_mid = (y_max + y_min) / 2.0
        # define origin point (where pattern starts)
        self.x_0 = x_mid + self.x_shift
        self.y_0 = y_mid + self.y_shift
        tape_paths, tape_angles, part_grid = tp.laminate_creation(
            tape_angles=self.t_angles, tape_widths=self.t_widths,
            tape_spacing=self.t_spacing, specimen_size=specimen_size,
            undulation_width=self.u_width, x_shift=self.x_shift,
            y_shift=self.y_shift)
        return tape_paths, tape_angles, part_grid


if __name__ == '__main__':
    # Material parameters
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 4 ply laminate

    # RVE dimensions
    x_min = -20.0
    x_max = 20.0
    y_min = -20.0
    y_max = 20.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = AP_PLY_Model('TestModel')
    mdl.time = 1.0  # set time attribute for testing purposes
    mdl.solid_mesh = 0.5  # set mesh size for testing purposes
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies, x_shift=10.0, y_shift=5.0)
    mdl.create_materials()
    paths_f, angles_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_job()
    mdl.save_model()
