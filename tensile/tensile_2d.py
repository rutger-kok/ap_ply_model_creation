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
from shapely.geometry import Polygon
from caeModules import *
from abaqus import *
from abaqusConstants import *
from ap_ply_elastic import AP_PLY_Elastic


class AP_PLY_Tensile_2D(AP_PLY_Elastic):
    def __init__(self, model_name, dir_name=None):
        AP_PLY_Elastic.__init__(self, model_name, dir_name)

    def set_test_parameters(self, time, test_speed, output_intervals,
                            shell_mesh, symmetry_mode):
        '''Set test parameters'''
        self.time = time
        self.test_speed = test_speed
        self.output_intervals = output_intervals
        self.shell_mesh = shell_mesh
        self.symmetry_mode = symmetry_mode

    def apply_constraints(self, symmetric=False):
        '''
        This method applies constraints to the assembly.

        Args:
            symmetric (boolean): False if the assembly does not contain
                mirrored instances (i.e. symmetry is enforced through a
                constraint rather than through geometric symmetry)

        Returns:
            None
        '''
        tol = 0.001
        x, y = zip(*self.shell_region.exterior.coords)
        x_max = max(x)  # determine max x-value for offsets
        y_max = max(y)
        x_min = min(x)  # determine max x-value for offsets
        y_min = min(y)
        z_min = 0.0
        z_max = self.l_thickness

        # identify edges at top and bottom for coupling
        top_edges_box = (x_min - tol, y_max - tol, z_min - tol, x_max + tol,
                         y_max + tol, z_max + tol)
        top_edges = self.get_edges_by_bounding_box(top_edges_box)
        top_edges_surface = self.assembly.Surface(side1Edges=top_edges,
                                                  name='Top Edges')
        bottom_edges_box = (x_min - tol, y_min - tol, z_min - tol, x_max + tol,
                            y_min + tol, z_max + tol)
        bottom_edges = self.get_edges_by_bounding_box(bottom_edges_box)
        bottom_edges_surface = self.assembly.Surface(side1Edges=bottom_edges,
                                                     name='Bottom Edges')
        symmetric_faces_box = (x_min - tol, y_min - tol, -tol, x_max + tol,
                               y_max + tol, tol)
        symmetric_faces = self.get_faces_by_bounding_box(symmetric_faces_box)
        symmetric_faces_set = self.assembly.Set(faces=symmetric_faces,
                                                name='Symmetry Faces')

        # create reference points for coupling
        rf_point_1_id = self.assembly.ReferencePoint(
            point=(0.0, y_max + 5.0, 0.0)).id
        rf_point_1 = self.assembly.referencePoints[rf_point_1_id]
        rf_point_1_region = self.assembly.Set(referencePoints=(rf_point_1,),
                                              name='Reference Point 1')
        rf_point_2_id = self.assembly.ReferencePoint(
            point=(0.0, y_min - 5.0, 0.0)).id
        rf_point_2 = self.assembly.referencePoints[rf_point_2_id]
        rf_point_2_region = self.assembly.Set(referencePoints=(rf_point_2,),
                                              name='Reference Point 2')

        self.model.Coupling(
            name='Top Coupling', controlPoint=rf_point_1_region,
            surface=top_edges_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

        self.model.Coupling(
            name='Bottom Coupling', controlPoint=rf_point_2_region,
            surface=bottom_edges_surface, influenceRadius=WHOLE_SURFACE,
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

        if symmetric is False:  # only create Z-symmetry if no mirrored parts
            self.model.ZsymmBC(name='Symmetry', createStepName='Loading Step',
                               region=symmetric_faces_set, localCsys=None)


if __name__ == '__main__':

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    output_intervals = 50  # requested field output intervals
    crosshead_velocity = 0.5
    symmetry_mode = 'BC'

    # Material parameters
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 8 ply laminate

    # Mesh sizes
    shell_mesh = 2.0

    # RVE dimensions
    x_min = -(tape_widths / 2.0)
    x_max = x_min + (tape_spacing + 1) * (tape_widths)
    y_min = -25.0
    y_max = 25.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = AP_PLY_Tensile_2D('TestModel')
    mdl.set_test_parameters(time, crosshead_velocity, output_intervals,
                            shell_mesh, symmetry_mode)
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies, x_shift=0.0, y_shift=0.0)
    mdl.create_materials()
    mdl.create_implicit_step()
    part = mdl.create_part_2d(specimen_size)
    grid, lines = mdl.create_grid_elastic(specimen_size)
    mdl.partition_part(part, lines, cells=False)
    mdl.place_tows_elastic(grid, specimen_size)
    mdl.assign_properties_2d(part, grid)
    mdl.mesh_part_2d(part)
    mdl.apply_constraints()
    mdl.create_job()
    mdl.save_model()
