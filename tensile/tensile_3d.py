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
from shapely.geometry import Polygon
from ap_ply_3d import AP_PLY_3D


class AP_PLY_Tensile_3D(AP_PLY_3D):
    def __init__(self, model_name, dir_name=None):
        AP_PLY_3D.__init__(self, model_name, dir_name)

    def set_test_parameters(self, time, test_speed, output_intervals,
                            solid_mesh, symmetry_mode):
        '''Set test parameters'''
        self.time = time  # simulation time in ms
        self.test_speed = test_speed  # crosshead velocity in mm/s
        self.output_intervals = output_intervals  # number of output intervals
        self.solid_mesh = solid_mesh  # mesh size in mm
        self.symmetry_mode = symmetry_mode  # 'BC' or 'GEOM'

    def apply_constraints(self):
        '''
        This method applies constraints to the assembly.

        Args:
            None

        Returns:
            None
        '''
        tol = 0.001
        x, y = zip(*self.specimen_size.exterior.coords)
        self.x_max = max(x)
        self.x_min = min(x)
        self.y_max = max(y)
        self.y_min = min(y)
        if self.symmetry_mode == 'BC':
            z_min = 0.0
            z_max = self.l_thickness / 2.0
        else:
            z_min = -self.l_thickness / 2.0
            z_max = self.l_thickness / 2.0

        # identify faces at top and bottom for coupling
        top_faces_box = (
            self.x_min - tol, self.y_max - tol, z_min - tol,
            self.x_max + tol, self.y_max + tol, z_max + tol)
        top_faces = self.get_faces_by_bounding_box(top_faces_box)
        top_faces_surface = self.assembly.Surface(
            side1Faces=top_faces, name='Top Faces')
        bottom_faces_box = (
            self.x_min - tol, self.y_min - tol, z_min - tol,
            self.x_max + tol, self.y_min + tol, z_max + tol)
        bottom_faces = self.get_faces_by_bounding_box(bottom_faces_box)
        bottom_faces_surface = self.assembly.Surface(
            side1Faces=bottom_faces, name='Bottom Faces')

        # create reference points for coupling
        rf_point_1_id = self.assembly.ReferencePoint(
            point=(0.0, self.y_max + 5.0, 0.0)).id
        rf_point_1 = self.assembly.referencePoints[rf_point_1_id]
        rf_point_1_region = self.assembly.Set(
            referencePoints=(rf_point_1,), name='Reference Point 1')
        rf_point_2_id = self.assembly.ReferencePoint(
            point=(0.0, self.y_min - 5.0, 0.0)).id
        rf_point_2 = self.assembly.referencePoints[rf_point_2_id]
        rf_point_2_region = self.assembly.Set(
            referencePoints=(rf_point_2,), name='Reference Point 2')

        self.model.Coupling(
            name='Top Coupling', controlPoint=rf_point_1_region,
            surface=top_faces_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

        self.model.Coupling(
            name='Bottom Coupling', controlPoint=rf_point_2_region,
            surface=bottom_faces_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
            u3=ON, ur1=ON, ur2=ON, ur3=ON)

        self.model.SmoothStepAmplitude(
            name='Smoothing Amplitude', timeSpan=STEP,
            data=((0.0, 0.0), (1e-05, 1.0)))

        self.model.DisplacementBC(
            name='Bottom Surface BC', createStepName='Loading Step',
            region=rf_point_2_region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)

        self.model.VelocityBC(
            name='Top Surface BC', createStepName='Loading Step',
            region=rf_point_1_region, v1=0.0, v2=self.test_speed, v3=0.0,
            vr1=0.0, vr2=0.0, vr3=0.0, amplitude='Smoothing Amplitude',
            localCsys=None, distributionType=UNIFORM, fieldName='')

        if self.symmetry_mode == 'BC':
            symmetric_faces_box = (
                self.x_min - tol, self.y_min - tol, -tol,
                self.x_max + tol, self.y_max + tol, tol)
            symmetric_faces = self.get_faces_by_bounding_box(
                symmetric_faces_box)
            symmetric_faces_set = self.assembly.Set(
                faces=symmetric_faces, name='Symmetry Faces')
            self.model.ZsymmBC(
                name='Symmetry', createStepName='Loading Step',
                region=symmetric_faces_set, localCsys=None)


if __name__ == '__main__':
    # Simulation parameters
    time = 2.5  # duration to simulate [ms]
    output_intervals = 100  # requested field output intervals
    crosshead_velocity = 0.5
    symmetry_mode = 'BC'

    # Material parameters
    tape_angles = (0, 45, 90, -45)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.205  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.1025  # ratio of undulation amplitude to length
    number_of_plies = 8  # symmetric 4 ply laminate

    # Mesh sizes
    solid_mesh = 1.0

    # RVE dimensions
    x_min = -20.0
    x_max = 20.0
    y_min = -20.0
    y_max = 20.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = AP_PLY_Tensile_3D('Test_QI')
    mdl.set_test_parameters(time, crosshead_velocity, output_intervals,
                            solid_mesh, symmetry_mode)
    mdl.create_explicit_step()
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies, x_shift=10.0, y_shift=10.0)
    mdl.create_materials()
    paths_f, angles_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_parts(grid_f)
    mdl.create_sequences()
    mdl.create_cohesive_interactions(tie=True)
    mdl.apply_constraints()
    mdl.create_job()
    mdl.save_model()
