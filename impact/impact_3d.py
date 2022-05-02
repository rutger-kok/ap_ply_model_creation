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
sys.path.append('C:\\Github\\ap_ply_model_creation\\impact')
from abaqus import *
from abaqusConstants import *
from shapely.geometry import Polygon
import mesh
import numpy as np
import regionToolset
from abaqus import *
from abaqusConstants import *
from shapely.geometry import Polygon
from impact_model import AP_PLY_Impact


class AP_PLY_Impact_Model_3D(AP_PLY_Impact):
    def __init__(self, model_name, dir_name=None):
        AP_PLY_Impact.__init__(self, model_name, dir_name)

    def create_3d_elastic_part(self, x_int, y_int, x_ext, y_ext):
        exterior = [(-x_ext, -y_ext), (x_ext, -y_ext), (x_ext, y_ext),
                    (-x_ext, y_ext), (-x_ext, -y_ext)]
        elastic_region = Polygon(exterior)

        if self.symmetry_mode == 'BC':
            z_min = 0.0
            z_max = self.l_thickness / 2.0
        else:
            z_min = -self.l_thickness / 2.0
            z_max = self.l_thickness / 2.0
        sketch = self.model.ConstrainedSketch(
            name='Elastic Part Sketch', sheetSize=200.0)
        sketch.rectangle(
            point1=(-x_ext, -y_ext), point2=(x_ext, y_ext))
        # Create a 3D deformable part
        part = self.model.Part(
            name='Elastic Solid', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseSolidExtrude(sketch=sketch, depth=(z_max - z_min))
        self.assembly.Instance(
            name='Elastic Solid Instance', part=part, dependent=ON)
        self.assembly.rotate(
            instanceList=('Elastic Solid Instance', ), angle=-90.0,
            axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0))
        self.assembly.translate(
            instanceList=('Elastic Solid Instance', ),
            vector=(0.0, -self.l_thickness / 2.0, 0.0))

        grid, lines = mdl.create_grid_elastic(elastic_region)
        mdl.partition_part(part, lines, cells=True)
        mdl.partition_through_thickness(part)
        mdl.place_tows_elastic(grid, elastic_region)
        mdl.assign_properties_3d(part, grid)

        # cut extrude interior out of part
        vertex_array = np.asarray([v.pointOn[0] for v in part.vertices])
        xmax, ymax, zmax = vertex_array.max(axis=0)
        xmin, ymin, zmin = vertex_array.min(axis=0)
        sketch_plane = part.faces.findAt(
            (xmax - 0.0001, ymin + 0.0001, zmax))
        sketch_edge = part.edges.findAt((xmax, ymin + 0.0001, zmax))
        t = part.MakeSketchTransform(
            sketchPlane=sketch_plane,
            sketchUpEdge=sketch_edge,
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
            origin=(self.x_0, self.y_0, zmax))
        s = self.model.ConstrainedSketch(
            name='Cutout-Elastic-1', sheetSize=200.0, gridSpacing=20.0,
            transform=t)
        s.CircleByCenterPerimeter(
            center=(self.x_0, self.y_0), point1=(x_int, 0.0))
        s.CircleByCenterPerimeter(
            center=(self.x_0, self.y_0), point1=(x_ext, 0.0))
        s.CircleByCenterPerimeter(
             center=(self.x_0, self.y_0), point1=(x_ext * 10, 0.0))
        part.CutExtrude(
            sketchPlane=sketch_plane, sketchUpEdge=sketch_edge,
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s,
            flipExtrudeDirection=OFF)
        self.mesh_part_3d(part)

    def create_shell_part(self, x_int, y_int, x_ext, y_ext):
        '''
        This function is used to create a (single!) 2D shell Tape part using
        geometric info from Shapely.
        '''
        exterior = [(x_ext, y_ext), (x_ext, -y_ext), (-x_ext, -y_ext),
                    (-x_ext, y_ext), (x_ext, y_ext)]
        shell_region = Polygon(exterior)
        x, y = zip(*shell_region.exterior.coords)

        # Create sketch
        sketch = self.model.ConstrainedSketch(
            name='2D Part Sketch', sheetSize=200.0)
        sketch.rectangle(point1=(-x_ext, -y_ext), point2=(x_ext, y_ext))
        part = self.model.Part(
            name='Shell Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseShell(sketch=sketch)

        instance_name = 'Shell Part Instance'
        self.assembly.Instance(name=instance_name, part=part, dependent=ON)
        self.assembly.rotate(
            instanceList=(instance_name, ), angle=-90.0,
            axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0))
        self.shell_edge_point = (x_ext, 0.0, 0.0)

        part_grid, lines = self.create_grid_elastic(shell_region)
        self.partition_part(part, lines, cells=False)
        self.place_tows_elastic(part_grid, shell_region)
        self.assign_properties_2d(part, part_grid)

        # cut sweep interior out of part
        hole_sketch = self.model.ConstrainedSketch(
            name='Cutout', sheetSize=200.0)
        vertex_array = np.asarray([v.pointOn[0] for v in part.vertices])
        xmax, ymax, zmax = vertex_array.max(axis=0)
        xmin, ymin, zmin = vertex_array.min(axis=0)
        sketch_plane = part.faces.findAt(
            (xmax - 0.0001, ymin + 0.0001, zmax))
        sketch_edge = part.edges.findAt((xmax, ymin + 0.0001, zmax))
        t = part.MakeSketchTransform(
            sketchPlane=sketch_plane,
            sketchUpEdge=sketch_edge,
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
            origin=(self.x_0, self.y_0, zmax))
        s = self.model.ConstrainedSketch(
            name='Sweep', sheetSize=200.0, gridSpacing=20.0, transform=t)
        s.retrieveSketch(sketch=self.model.sketches['Cutout'])
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(x_int, 0.0))
        part.CutExtrude(
            sketchPlane=sketch_plane, sketchUpEdge=sketch_edge,
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s,
            flipExtrudeDirection=OFF)
        self.mesh_part_2d(part)

    def trim_solid_circular(self, x_int, y_int):
        # trim solid part to make it circular
        for i, part in enumerate(self.model.parts.values()):
            if 'Layer' in part.name:
                vertex_array = np.asarray(
                    [v.pointOn[0] for v in part.vertices])
                xmax, ymax, zmax = vertex_array.max(axis=0)
                xmin, ymin, zmin = vertex_array.min(axis=0)
                sketch_plane = part.faces.findAt(
                    (xmax - 0.0001, ymin + 0.0001, zmax))
                sketch_edge = part.edges.findAt((xmax, ymin + 0.0001, zmax))
                t = part.MakeSketchTransform(
                    sketchPlane=sketch_plane,
                    sketchUpEdge=sketch_edge,
                    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
                    origin=(0.0, 0.0, zmax))
                s = self.model.ConstrainedSketch(
                    name='Cutout-{}'.format(i), sheetSize=200.0,
                    gridSpacing=20.0, transform=t)
                s.CircleByCenterPerimeter(
                    center=(self.x_0, self.y_0), point1=(x_int, 0.0))
                s.rectangle(
                    point1=(-500.0, -500.0), point2=(500.0, 500.0))
                part.CutExtrude(
                    sketchPlane=sketch_plane, sketchUpEdge=sketch_edge,
                    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s,
                    flipExtrudeDirection=OFF)
                self.mesh_part_3d(part)

    def rotate_3d_instances(self):
        non_3d_instances = ['Impactor Instance', 'Bottom Plate Instance',
                            'Clamp Instance 1', 'Clamp Instance 2',
                            'Clamp Instance 3', 'Clamp Instance 4',
                            'Shell Part Instance', 'Elastic Solid Instance']
        instances_3d = [inst.name for inst in self.assembly.instances.values()
                        if inst.name not in non_3d_instances]
        self.assembly.rotate(instanceList=instances_3d, angle=-90.0,
                             axisPoint=(0.0, 0.0, 0.0),
                             axisDirection=(1.0, 0.0, 0.0))

    def elastic_region_tie(self, x_int, y_int):
        # determine surfaces to couple
        damage_instances = [
            inst for inst in self.assembly.instances.values()
            if 'Layer' in inst.name]
        # define tolerance for search
        tol = 1.0e-3
        i = 0
        # search by checking if point on face/edge is on circle
        for inst in damage_instances:
            for face in inst.faces:
                x, y, z = face.pointOn[0]
                eq = x**2.0 + z**2.0  # z because part is rotated in assembly
                if abs(eq - x_int**2.0) <= tol:
                    if i == 0:
                        tie_faces = inst.faces.findAt(
                            (face.pointOn[0], ))
                    else:
                        tie_faces = tie_faces + inst.faces.findAt(
                            (face.pointOn[0], ))
                    i += 1

        # determine shell edges to couple
        elastic_instance = self.assembly.instances['Elastic Solid Instance']
        j = 0
        for face in elastic_instance.faces:
            x, y, z = face.pointOn[0]
            eq = x**2.0 + z**2.0
            if abs(eq - x_int**2.0) <= tol:
                if j == 0:
                    tie_faces_elastic = elastic_instance.faces.findAt(
                        (face.pointOn[0], ))
                else:
                    tie_faces_elastic = tie_faces_elastic + elastic_instance.faces.findAt(
                        (face.pointOn[0], ))
                j += 1

        # define surfaces for coupling
        tie_region_elastic = self.assembly.Surface(
            side1Faces=tie_faces_elastic, name='Tie Faces Elastic')        
        tie_region_damage = self.assembly.Surface(
            side1Faces=tie_faces, name='Tie Faces Damage')
        self.model.Tie(
            name='Elastic Tie', master=tie_region_damage, adjust=ON,
            slave=tie_region_elastic, positionToleranceMethod=COMPUTED,
            tieRotations=ON, thickness=ON)

    def solid_shell_coupling(self, x_int, y_int):
        # define tolerance for search
        tol = 1.0e-3
        i = 0
        # search by checking if point on face/edge is on circle
        elastic_instance = self.assembly.instances['Elastic Solid Instance']
        for face in elastic_instance.faces:
            x, y, z = face.pointOn[0]
            eq = x**2.0 + z**2.0  # z because part is rotated in assembly
            if abs(eq - x_int**2.0) <= tol:
                if i == 0:
                    couple_faces = elastic_instance.faces.findAt(
                        (face.pointOn[0], ))
                else:
                    couple_faces = couple_faces + elastic_instance.faces.findAt(
                        (face.pointOn[0], ))
                i += 1

        # determine shell edges to couple
        shell_instance = self.assembly.instances['Shell Part Instance']
        j = 0
        for edge in shell_instance.edges:
            x, y, z = edge.pointOn[0]
            eq = x**2.0 + z**2.0
            if abs(eq - x_int**2.0) <= tol:
                if j == 0:
                    couple_edges = shell_instance.edges.findAt(
                        (edge.pointOn[0], ))
                else:
                    couple_edges = couple_edges + shell_instance.edges.findAt(
                        (edge.pointOn[0], ))
                j += 1

        # define surfaces for coupling
        coupling_region_shell = self.assembly.Surface(
            side1Edges=couple_edges, name='Coupling Edges Shell')        
        coupling_region_solid = self.assembly.Surface(
            side1Faces=couple_faces, name='Coupling Faces Solid')
        self.model.ShellSolidCoupling(
            name='Shell-Solid-Coupling', shellEdge=coupling_region_shell,
            solidFace=coupling_region_solid, positionToleranceMethod=COMPUTED)


    def create_cohesive_interactions(self, tie=False):
        # create explicit general contact definition
        self.model.ContactExp(name='GC', createStepName='Loading Step')
        gc = self.model.interactions['GC']
        gc.includedPairs.setValuesInStep(
            stepName='Loading Step', useAllstar=ON)
        # define 'Tangential' as default behaviour
        gc.contactPropertyAssignments.appendInStep(
            stepName='Loading Step',
            assignments=((GLOBAL, SELF, 'Tangential'), ))
        # determine contacts\
        non_cohesive_instances = [
            'Impactor Instance', 'Bottom Plate Instance',
            'Clamp Instance 1', 'Clamp Instance 2', 'Clamp Instance 3',
            'Clamp Instance 4', 'Shell Part Instance',
            'Elastic Solid Instance']
        cohesive_instances = []
        for inst in self.assembly.instances.values():
            if inst.name not in non_cohesive_instances:
                cohesive_instances.append(inst.name)
        if tie:
            self.model.contactDetection(
                defaultType=TIE, nameEachSurfaceFound=OFF,
                searchDomain=cohesive_instances, createUnionOfSlaveSurfaces=ON,
                createUnionOfMasterSurfaces=ON, separationTolerance=0.0001)
        else:
            self.model.contactDetection(
                defaultType=CONTACT, interactionProperty='Cohesive',
                nameEachSurfaceFound=OFF, createUnionOfMasterSurfaces=ON,
                createUnionOfSlaveSurfaces=ON, searchDomain=cohesive_instances,
                separationTolerance=0.0001)
            # assign cohesive behaviour to contacting tape surfaces
            for inter in self.model.interactions.values():
                if inter.name == 'GC':
                    continue
                else:
                    master_name = inter.master[0]
                    slave_name = inter.slave[0]
                    master_surface = self.assembly.surfaces[master_name]
                    slave_surface = self.assembly.surfaces[slave_name]
                    gc.contactPropertyAssignments.appendInStep(
                        stepName='Loading Step',
                        assignments=(
                            (master_surface, slave_surface, 'Cohesive'), ))
                    del self.model.interactions['{}'.format(inter.name)]
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
        clamp_instance_keys = [
            'Clamp Instance 1', 'Clamp Instance 2', 'Clamp Instance 3',
            'Clamp Instance 4']
        for i, clamp_instance in enumerate(clamp_instance_keys):
            instance_faces = self.assembly.instances[clamp_instance].faces
            bottom_face = instance_faces.getByBoundingBox(
                -200.0, self.l_thickness / 2.0 - 0.01, -200.0,
                200.0, self.l_thickness / 2.0 + 0.01, 200.0)
            clamp_face = self.assembly.Surface(
                side1Faces=bottom_face, name='Clamp Face {}'.format(i + 1))
            self.model.SurfaceToSurfaceContactExp(
                name='Shell-Clamp-{}'.format(i + 1), master=clamp_face,
                createStepName='Loading Step', slave=shell_top_face,
                mechanicalConstraint=KINEMATIC, sliding=FINITE,
                interactionProperty='Tangential', initialClearance=OMIT,
                datumAxis=None, clearanceRegion=None)
            gc.excludedPairs.setValuesInStep(
                stepName='Loading Step',
                addPairs=((clamp_face, shell_top_face), ))
        self.model.SurfaceToSurfaceContactExp(
            name='Shell-Support', master=shell_top_face,
            createStepName='Loading Step', slave=support_top_face,
            mechanicalConstraint=KINEMATIC, sliding=FINITE,
            interactionProperty='Tangential', initialClearance=OMIT,
            datumAxis=None, clearanceRegion=None)

        # exclude contact between surfaces at tie and coupling constraints
        surf1 = self.assembly.surfaces['Tie Faces Elastic']
        surf2 = self.assembly.surfaces['Tie Faces Damage']
        surf3 = self.assembly.surfaces['Coupling Edges Shell']
        surf4 = self.assembly.surfaces['Coupling Faces Solid']
        gc.excludedPairs.setValuesInStep(
            stepName='Loading Step',
            addPairs=((surf1, surf2), (surf3, surf4),
                      (shell_top_face, support_top_face)))


if __name__ == '__main__':

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    output_intervals = 50  # requested field output intervals
    energy = 30.0  # impact energy to simulate
    symmetry_mode = 'GEOM'

    # Material parameters
    # tape_angles = (0, 90)  # define angles of tapes
    # tape_angles = (0, 45, -45, 90)  # define angles of tapes
    tape_angles = (0, 60, -60)  # define angles of tapes
    tape_widths = 12.7
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.213  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.104  # ratio of undulation amplitude to length
    number_of_plies = 24  # symmetric 8 ply laminate

    # Mesh sizes
    solid_mesh = (cured_ply_thickness / undulation_ratio) / 2.0
    shell_mesh = 1.5

    # Define specimen dimensions
    # Fine mesh region
    xMin_f = -31.0
    xMax_f = -xMin_f
    yMin_f = -31.0
    yMax_f = -yMin_f
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
    # Create model
    mdl = AP_PLY_Impact_Model_3D('DWT_60_baseline_full_test')
    mdl.set_test_parameters(
        time, output_intervals, energy, shell_mesh, solid_mesh, symmetry_mode)
    mdl.set_specimen_parameters(
        tape_angles, tape_widths, tape_spacing, cured_ply_thickness,
        undulation_ratio, number_of_plies, x_shift=0.0, y_shift=0.0)
    mdl.create_materials(material_name='HiTapeUD210')
    # create experimental apparatus parts
    mdl.create_impactor_part()
    mdl.create_clamp_part()
    mdl.create_bottom_plate()
    # create 3D solid parts
    paths_f, angles_f, grid_f = mdl.create_tape_paths(fine_region)
    mdl.create_parts(grid_f)
    mdl.create_sequences()
    mdl.trim_solid_circular(x_int_el, y_int_el)
    # create shell part
    mdl.create_shell_part(x_int, y_int, x_ext, y_ext)
    # create elastic part
    mdl.create_3d_elastic_part(x_int_el, y_int_el, x_ext_el, y_ext_el)
    # # create step
    mdl.create_explicit_step()
    # # populate assembly with parts and rotate to correct orientation
    mdl.create_assembly()
    mdl.rotate_3d_instances()
    # # apply constraints
    mdl.constrain_b_plate()
    mdl.constrain_clamps()
    mdl.constrain_impactor()
    # # tie shell parts to 3D solid parts
    mdl.solid_shell_coupling(x_int, y_int)
    mdl.elastic_region_tie(x_int_el, y_int_el)
    # # create tie constraints or cohesive interactions between layers
    mdl.create_cohesive_interactions()
    # # create job and save model
    mdl.create_job()
    mdl.save_model()
