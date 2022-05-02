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
sys.path.append('C:\\Github\\interlaced_model_creation\\impact')
from abaqus import *
from abaqusConstants import *
from shapely.geometry import Polygon
import mesh
import regionToolset
from abaqus import *
from abaqusConstants import *
from shapely.geometry import Polygon
from impact_model_quarter import AP_PLY_Impact_quarter


class AP_PLY_Impact_Model_2D(AP_PLY_Impact_quarter):
    def __init__(self, model_name, dir_name=None):
        AP_PLY_Impact_quarter.__init__(self, model_name, dir_name)

    def create_shell_part(self, x_int, y_int, x_ext, y_ext):
        '''
        This method is used to create a (single!) 2D shell Tape part using
        geometric info from Shapely.
        '''
        exterior = [(0.0, 0.0), (x_ext, 0.0), (x_ext, y_ext),
                    (0.0, y_ext), (0.0, 0.0)]
        shell_region = Polygon(exterior)

        # define origin point (where pattern starts)
        self.x_0 = 0.0
        self.y_0 = 0.0

        # Create sketch
        sketch = self.model.ConstrainedSketch(
            name='2D Part Sketch', sheetSize=200.0)
        sketch.rectangle(point1=(0.0, 0.0), point2=(x_ext, y_ext))
        part = self.model.Part(
            name='Shell Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseShell(sketch=sketch)

        instance_name = 'Shell Part Instance'
        self.assembly.Instance(name=instance_name, part=part, dependent=ON)

        part_grid, lines = self.create_grid_elastic(shell_region)
        self.partition_part(part, lines, cells=False)
        self.place_tows_elastic(part_grid, shell_region)
        self.assign_properties_2d(part, part_grid)
        self.mesh_part_2d(part)

    def rotate_3d_instances(self):
        non_3d_instances = ['Impactor Instance', 'Bottom Plate Instance',
                            'Clamp Instance']
        instances_3d = [inst.name for inst in self.assembly.instances.values()
                        if inst.name not in non_3d_instances]
        self.assembly.rotate(instanceList=instances_3d, angle=-90.0,
                             axisPoint=(0.0, 0.0, 0.0),
                             axisDirection=(1.0, 0.0, 0.0))

    def apply_symmetry(self):
        tol = 0.0001
        shell_instance = self.assembly.instances['Shell Part Instance']
        # find faces/edges on part instances for symmetry
        # z-symmetry
        shell_z_sym = shell_instance.edges.getByBoundingBox(
            -1.0, -self.l_thickness - tol, 0.0 - tol,
            75.0, self.l_thickness + tol, 0.0 + tol)
        z_sym_region = self.assembly.Set(
            edges=shell_z_sym, name='Z-Symmetry')
        self.model.ZsymmBC(
            name='Z-Symmetry', createStepName='Initial',
            region=z_sym_region, localCsys=None)
        # x-symmetry
        shell_x_sym = shell_instance.edges.getByBoundingBox(
            0.0 - tol, -self.l_thickness - tol, -50.0,
            0.0 + tol, self.l_thickness + tol, 0.0 + tol)
        x_sym_region = self.assembly.Set(
            edges=shell_x_sym, name='X-Symmetry')
        self.model.XsymmBC(
            name='X-Symmetry', createStepName='Initial', region=x_sym_region,
            localCsys=None)

    def create_cohesive_interactions(self, tie=False):
        # create contact definition
        # determine shell edges to couple
        impactor_instance = self.assembly.instances['Impactor Instance']
        j = 0
        tol = 1.0e-3
        for face in impactor_instance.faces:
            x, y, z = face.pointOn[0]
            eq1 = x**2.0 + (y - 10.566)**2.0 + z**2.0
            eq2 = x**2.0 + z**2.0
            if abs(eq1 - 8.0**2.0) <= tol or (abs(eq2 - 8.0**2.0) <= tol):
                if j == 0:
                    impactor_faces = impactor_instance.faces.findAt(
                        (face.pointOn[0], ))
                else:
                    impactor_faces = impactor_faces + impactor_instance.faces.findAt(
                        (face.pointOn[0], ))
                j += 1
        impactor_surface = self.assembly.Surface(
            side1Faces=impactor_faces, name='Impactor Faces')
        # define contact between support and shell region, as well as shell
        # region and clamp, separately
        sif = self.get_faces_by_bounding_box(
            [-200.0, -0.01, -200.0, 200.0, 0.01, 200.0],
            instances=[self.assembly.instances['Shell Part Instance'], ])  
        shell_top_face = self.assembly.Surface(
            side1Faces=sif,
            name='Shell Faces')
        bif = self.assembly.instances['Bottom Plate Instance'].faces
        support_top_face = self.assembly.Surface(
            side1Faces=bif.findAt(((22.5, -self.l_thickness / 2.0, -50.0), )),
            name='Support Faces')
        cif = self.assembly.instances['Clamp Instance'].faces
        clamp_bottom_face = self.assembly.Surface(
            side1Faces=cif.findAt(((50.0, self.l_thickness / 2.0, -44.0), )),
            name='Clamp Bottom Face')

        self.model.SurfaceToSurfaceContactExp(
            name='Shell-Clamp', master=clamp_bottom_face,
            createStepName='Loading Step', slave=shell_top_face,
            mechanicalConstraint=KINEMATIC, sliding=FINITE,
            interactionProperty='Tangential', initialClearance=OMIT,
            datumAxis=None, clearanceRegion=None)
        self.model.SurfaceToSurfaceContactExp(
            name='Impactor-Shell', master=impactor_surface,
            createStepName='Loading Step', slave=shell_top_face,
            mechanicalConstraint=KINEMATIC, sliding=FINITE,
            interactionProperty='Tangential', initialClearance=OMIT,
            datumAxis=None, clearanceRegion=None)
        self.model.SurfaceToSurfaceContactExp(
            name='Shell-Support', master=shell_top_face,
            createStepName='Loading Step', slave=support_top_face,
            mechanicalConstraint=KINEMATIC, sliding=FINITE,
            interactionProperty='Tangential', initialClearance=OMIT,
            datumAxis=None, clearanceRegion=None)


if __name__ == '__main__':

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    output_intervals = 50  # requested field output intervals
    energy = 30.0  # impact energy to simulate
    symmetry_mode = 'GEOM'

    # Material parameters
    tape_angles = (0, 45, -45, 90)  # define angles of tapes
    # tape_angles = (0, 90)  # define angles of tapes
    # tape_angles = (0, 60, -60)  # define angles of tapes
    tape_widths = 12.7
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.213  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.104  # ratio of undulation amplitude to length
    number_of_plies = 24  # symmetric 8 ply laminate

    # Mesh sizes
    solid_mesh = 1.5
    shell_mesh = 1.5

    # Define specimen dimensions
    # Fine mesh region
    xMin_f = 0.0
    xMax_f = 29.0
    yMin_f = 0.0
    yMax_f = 29.0
    fine_region = Polygon([(xMin_f, yMin_f), (xMax_f, yMin_f),
                           (xMax_f, yMax_f), (xMin_f, yMax_f)])

    # Shell region
    x_int = 36.5
    x_ext = 75.0
    y_int = 36.5
    y_ext = 50.0

    # Elastic region
    x_int_el = xMax_f
    x_ext_el = 36.5
    y_int_el = yMax_f
    y_ext_el = 36.5

    shift_x = -xMax_f / 2.0
    shift_y = -yMax_f / 2.0
    # Create model
    mdl = AP_PLY_Impact_Model_2D('DWT_QI_quarter_2d')
    mdl.set_test_parameters(
        time, output_intervals, energy, shell_mesh, solid_mesh, symmetry_mode)
    mdl.set_specimen_parameters(
        tape_angles, tape_widths, tape_spacing, cured_ply_thickness,
        undulation_ratio, number_of_plies, x_shift=shift_x, y_shift=shift_y)
    mdl.create_materials(material_name='HiTapeUD210')
    # create experimental apparatus parts
    mdl.create_impactor_part()
    mdl.create_clamp_part()
    mdl.create_bottom_plate()
    # create shell part
    mdl.create_shell_part(x_int, y_int, x_ext, y_ext)
    # create step
    mdl.create_explicit_step()
    # # populate assembly with parts and rotate to correct orientation
    mdl.create_assembly()
    mdl.rotate_3d_instances()
    # # apply constraints
    mdl.constrain_b_plate()
    mdl.constrain_clamps()
    mdl.constrain_impactor()
    mdl.apply_symmetry()
    mdl.create_cohesive_interactions()
    # # create job and save model
    mdl.create_job()
    mdl.save_model()
