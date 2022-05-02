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
from caeModules import *
import mesh
import regionToolset
from shapely.geometry import Polygon
from shapely.ops import unary_union
from ap_ply_model import AP_PLY_Model
import numpy as np
from random import shuffle
from contextlib import contextmanager


@contextmanager
def stdout_redirected(new_stdout):
    '''
    Function used to redicrect the standard error output to a new location.
    '''
    save_stdout = sys.stdout
    sys.stdout = new_stdout
    try:
        yield None
    finally:
        sys.stdout = save_stdout


class AP_PLY_3D(AP_PLY_Model):
    def __init__(self, model_name, dir_name=None):
        AP_PLY_Model.__init__(self, model_name, dir_name)

    def group_by_type(self, part_num, grid, n, obj_type, angle, und_pair=None):
        '''
        Method used to group the objects in a grid of shapely polygons by their
        type, angle, and (if they are undulation regions) the angles of the
        undulation constituents.
        '''
        if und_pair is None:
            same_type = [
                obj for obj in grid.values()
                if obj.object_type[n - 1] == obj_type
                and obj.angle[n - 1] == angle]
        else:
            same_type = [
                obj for obj in grid.values()
                if obj.object_type[n - 1] == obj_type
                and obj.angle[n - 1] == angle
                and obj.und_pairs[n - 1] == und_pair]

        merged_obj = unary_union(same_type)
        # split merged_obj into  its constituent Polygons
        polygons = self.split_multipolygon(merged_obj)

        # iterate over polygons to create parts
        parts = []
        for poly in polygons:
            if isinstance(und_pair, int):
                part_id = '{}_{}_{}_{}_{}'.format(
                    obj_type, angle, n, und_pair, part_num)
            elif isinstance(und_pair, tuple):
                part_id = '{}_{}_{}_{}_{}'.format(
                    obj_type, angle, n, '-'.join(str(x) for x in und_pair),
                    part_num)
            else:
                part_id = '{}_{}_{}_{}'.format(
                    obj_type, angle, n, part_num)
            part, orient_angle = self.create_3d_part(
                poly, part_id, obj_type, n, angle, und_pair)
            parts.append((
                part_id, part, orient_angle,
                poly.representative_point().coords[0]))
            part_num += 1

        instances = []
        # create part instances
        for part_id, part, orient_angle, centroid in parts:
            instance = self.create_instances(
                part, part_id, orient_angle, centroid, n)
            instances.append(instance)
        return instances, part_num

    def create_parts(self, grid, mesh=False):
        '''
        This function is used to create parts using geometric info from
        Shapely.

        Args:
            grid (dict): dictionary containing Shapely polygons

        Returns:
            None
        '''
        region_types = ['Tape', 'Resin', 'Undulation']
        und_pairs = set()
        for obj in grid.values():
            for i in range(1, len(self.t_angles) + 1):
                if obj.object_type[i - 1] == 'Undulation':
                    und_pairs.add(obj.und_pairs[i - 1])
        und_pairs = list(und_pairs)
        part_num = 0
        for n in range(1, len(self.t_angles) + 1):
            layer_instances = []  # instances in layer
            name = 'Layer-{}'.format(n)  # layer name
            # define dictionary to store material orientations
            self.mat_orient = {name: {}}
            for obj_type in region_types:
                for angle in self.t_angles:
                    if obj_type == 'Undulation':
                        for und_pair in und_pairs:
                            inst, part_num = self.group_by_type(
                                part_num, grid, n, obj_type, angle,
                                und_pair=und_pair)
                            layer_instances.extend(inst)
                    else:
                        und_pair = None
                        inst, part_num = self.group_by_type(
                            part_num, grid, n, obj_type, angle)
                        layer_instances.extend(inst)
            # merge part instances in a single layer into one part
            self.merge_layer(name, layer_instances, mesh=mesh)

    def shuffle_merge(self, layer_name, instances):
        '''
        Abaqus does not reaise warnings according to the standard Python
        protocol. As such, to catch warnings related to incomplete merging
        of part instances, redirect stdout to a file temporarily and check
        that file for warnings. If a warning is present, shuffle instances
        and retry merging. This process is recursive.
        '''
        warning_flag = 0
        warning_file = 'warnings.txt'
        with open(warning_file, "w+") as f:
            with stdout_redirected(f):
                shuffle(instances)
                self.assembly.InstanceFromBooleanMerge(
                    name=layer_name, keepIntersections=ON, instances=instances,
                    domain=GEOMETRY, originalInstances=DELETE)
                content = f.readlines()
                if 'Warning' in content:
                    warning_flag = 1
        if warning_flag:
            self.shuffle_merge(layer_name, instances)

    def merge_layer(self, layer_name, instances, mesh=False):
        '''
        Merge objects in the same layer. When merging fails shuffling the order
        of the instances to be merged often fixes the issue.
        '''
        delete_parts_list = [inst.partName for inst in instances]
        self.shuffle_merge(layer_name, instances)
        merged_part = self.model.parts[layer_name]
        self.assign_orientations(merged_part, layer_name)
        if mesh:
            self.mesh_part_3d(merged_part)
        for part_name in delete_parts_list:
            del self.model.parts[part_name]

    def assign_orientations(self, merged_part, layer_name):
        '''
        Assign material orientations to the merged parts.
        '''
        cells = merged_part.cells
        for theta, regions in self.mat_orient[layer_name].iteritems():
            for i, rep_point in enumerate(regions):
                if i == 0:
                    selected_cells = cells.findAt((rep_point, ))
                else:
                    selected_cells += cells.findAt((rep_point, ))
            # assign orientation to region
            region = merged_part.Set(
                cells=selected_cells, name='Cells-{}'.format(int(theta)))
            merged_part.MaterialOrientation(
                region=region, orientationType=SYSTEM, localCsys=None,
                stackDirection=STACK_3, additionalRotationField='',
                angle=theta, additionalRotationType=ROTATION_ANGLE,
                axis=AXIS_3, fieldName='')

    def create_sequences(self):
        '''
        If the AP-PLY laminate consists of multiple "sets" of AP-PLY
        sublaminates, then this method creates the required instances.
        '''
        seq_offset = len(self.t_angles) * self.t_thickness
        part_names = [
            p for p in self.model.parts.keys() if 'Layer' in p]
        for sequence in range(self.sequences - 1):
            # create new instances:
            for name in part_names:
                if 'Layer' in name:
                    part = self.model.parts[name]
                    new_instance = self.assembly.Instance(
                        name=name + '-{}'.format(sequence + 2), part=part,
                        dependent=ON)
                    self.assembly.translate(
                        instanceList=(new_instance.name, ),
                        vector=(0.0, 0.0, (sequence + 1) * seq_offset))
        if self.symmetry_mode == 'GEOM':
            for inst in self.assembly.instances.values():
                if 'Layer' in inst.name:
                    part = self.model.parts[inst.partName]
                    sym_instance = self.assembly.Instance(
                        name=inst.name + '-sym', part=part, dependent=ON)
                    layer_number = int(
                        inst.name[inst.name.find('-') + 1: inst.name.rfind('-')])
                    sequence = int(inst.name[inst.name.rfind('-') + 1:])
                    layer_z = (2 * layer_number - 1) * self.t_thickness
                    seq_z = (sequence - 1) * seq_offset
                    self.assembly.translate(
                        instanceList=(sym_instance.name, ),
                        vector=(0.0, 0.0, -1.0 * (layer_z + seq_z)))

    def create_cohesive_interactions(self, tie=False):
        '''
        Create cohesive interactions between the plies in an AP-PLY laminate.
        If "Tie" is specified then the plies will be connected to eachother
        using a Tie constraint, otherwise cohesive surfaces will be defined.

        Note: In an Abaqus Explicit simulation, cohesive contacts must be
        defined as part of a general contact definition, they cannot be defined
        as surface-to-surface contacts. That is why in this method the
        contactDetection built-in is used to find the surfaces that are in
        contact, and then these interactions are converted to individual
        contacts in the general contact definition.
        '''
        if tie:
            self.model.contactDetection(
                defaultType=TIE, nameEachSurfaceFound=OFF, searchDomain=MODEL,
                createUnionOfMasterSurfaces=ON, separationTolerance=0.0001,
                createUnionOfSlaveSurfaces=ON)
        else:
            self.model.contactDetection(
                defaultType=CONTACT, interactionProperty='Cohesive',
                nameEachSurfaceFound=OFF, createUnionOfMasterSurfaces=ON,
                createUnionOfSlaveSurfaces=ON, searchDomain=MODEL,
                separationTolerance=0.0001)
            # create explicit general contact definition
            self.model.ContactExp(name='GC', createStepName='Initial')
            gc = self.model.interactions['GC']
            gc.includedPairs.setValuesInStep(
                stepName='Initial', useAllstar=ON)
            # define 'Tangential' as default behaviour
            gc.contactPropertyAssignments.appendInStep(
                stepName='Initial',
                assignments=((GLOBAL, SELF, 'Tangential'), ))
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
                        stepName='Initial',
                        assignments=(
                            (master_surface, slave_surface, 'Cohesive'), ))
                    del self.model.interactions['{}'.format(inter.name)]

    def split_multipolygon(self, geom):
        '''
        Recursively split a MultiPolygon or GeometryCollection into individual
        Polygons.
        '''
        if geom.geom_type == 'Polygon':
            yield geom
        else:
            for sub_geom in geom:
                if (sub_geom.geom_type == 'GeometryCollection' or
                    sub_geom.geom_type == 'MultiPolygon'):
                        for sub_sub_geom in sub_geom:
                            yield self.split_multipolygon(sub_sub_geom)
                elif sub_geom.geom_type == 'Polygon':
                    yield sub_geom
                else:
                    print 'Error:', sub_geom.geom_type

    def sketch_part(self, obj, obj_id):
        '''
        Method to create a sketch of a part using the coordinates from
        Shapely.

        Args:
            obj (Shapely Polygon): the polygon to be sketched.
            obj_id (int): the numeric ID of the polygon.

        Returns:
            sketch (Abaqus sketch): sketch of the obj polygon.
        '''
        sketch = self.model.ConstrainedSketch(
            name='Sketch {}'.format(obj_id), sheetSize=200.0)
        # plot interior points of object
        for hole in obj.interiors:
            x_int, y_int = hole.xy
            abaqus_coords_int = zip(x_int, y_int)
            for idx in range(len(abaqus_coords_int) - 1):
                sketch.Line(point1=(abaqus_coords_int[idx][0],
                                    abaqus_coords_int[idx][1]),
                            point2=(abaqus_coords_int[idx + 1][0],
                                    abaqus_coords_int[idx + 1][1]))

        # clean Polygon vertices and plot exterior points
        clean_obj = self.simplify_polygon(obj)
        x_ext, y_ext = clean_obj.exterior.xy
        abaqus_coords_ext = zip(x_ext, y_ext)
        for idx in range(len(abaqus_coords_ext) - 1):
            sketch.Line(point1=(abaqus_coords_ext[idx][0],
                                abaqus_coords_ext[idx][1]),
                        point2=(abaqus_coords_ext[idx + 1][0],
                                abaqus_coords_ext[idx + 1][1]))
        return sketch

    def create_3d_part(self, obj, obj_id, obj_type, obj_layer, obj_angle,
                       und_pairs):
        '''
        Method to create a 3D part instance and assign properties to it.
        '''
        part = self.model.Part(
            name='{}'.format(obj_id), dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        try:
            sketch = self.sketch_part(obj, obj_id)
            part.BaseSolidExtrude(sketch=sketch, depth=self.t_thickness)
        except:  # catch self-intersect error
            eps = 1.0e-8
            obj = obj.buffer(eps).buffer(-eps)
            sketch = self.sketch_part(obj, obj_id)
            try:
                part.BaseSolidExtrude(sketch=sketch, depth=self.t_thickness)
            except:
                print obj
                raise ValueError
        cell_region = regionToolset.Region(cells=part.cells[:])
        if obj_type == 'Tape':
            section_name = 'Tape-Damage'
            orient_angle = obj_angle
        elif obj_type == 'Resin':
            section_name = 'Resin-Rich-Damage'
            orient_angle = obj_angle
        elif obj_type == 'Undulation':
            undulation_angles = und_pairs
            if isinstance(undulation_angles, int):
                section_name = "Undulation-Damage-(0, 'Resin')"
                orient_angle = obj_angle
            else:
                section_angle = (undulation_angles[0], undulation_angles[1])
                section_name = 'Undulation-Damage-{}'.format(section_angle)
                orient_angle = undulation_angles[1]
        part.SectionAssignment(
            region=cell_region, sectionName=section_name,
            offsetType=MIDDLE_SURFACE, thicknessAssignment=FROM_SECTION,
            offset=0.0, offsetField='')
        return part, orient_angle

    def create_instances(self, part, obj_id, obj_angle, rep_point, obj_layer):
        '''
        Create part instances (from merged layer) and offset them by the
        required value in the z-drection.
        '''
        # create the part instance
        z_offset = (obj_layer - 1) * self.t_thickness
        instance = self.assembly.Instance(
            name=obj_id, part=part, dependent=ON)
        self.assembly.translate(
            instanceList=(obj_id, ), vector=(0.0, 0.0, z_offset))
        # define material orientation
        z_coord = self.t_thickness * (2 * obj_layer - 1) / 2.0
        self.mat_orient['Layer-{}'.format(obj_layer)].setdefault(
            obj_angle, []).append((rep_point[0], rep_point[1], z_coord))
        return instance

    def get_angles(self, vec_1, vec_2):
        '''
        Returns the angle, in degrees, between two vectors
        Adapted from: https://github.com/Toblerity/Shapely/issues/1046
        '''
        dot = np.dot(vec_1, vec_2)
        det = np.cross(vec_1, vec_2)
        angle_in_rad = np.arctan2(det, dot)
        return np.degrees(angle_in_rad)

    def simplify_polygon(self, poly_in, deg_tol=1.0):
        '''
        Remove "unnecessary" vertices from Polygon by calculating the change
        in direction between each vertex and removing vertices for which the
        change in direction is 0.

        Args:
            poly_in (Shapely Polygon): the polygon to be simplified
            deg_tol (float): degree tolerance for comparison between
                successive vectors

        Returns:
            clean_polygon (Shapely Polygon): 'simplified' polygon

        Adapted from: https://github.com/Toblerity/Shapely/issues/1046
        '''
        ext_poly_coords = poly_in.exterior.coords[:]
        vector_rep = np.diff(ext_poly_coords, axis=0)
        angles_list = []
        for i in range(0, len(vector_rep) - 1):
            angles_list.append(
                np.abs(self.get_angles(vector_rep[i], vector_rep[i + 1])))
        # add angle between last and first vector
        angles_list.append(
            np.abs(self.get_angles(vector_rep[0], vector_rep[-1])))
        # get mask satisfying tolerance
        thresh_vals_by_deg = (
            np.where(np.array(angles_list) > deg_tol)[0] + 1).tolist()
        new_vertices = [ext_poly_coords[idx] for idx in thresh_vals_by_deg]
        clean_polygon = Polygon(new_vertices)
        return clean_polygon

    def mesh_part_3d(self, part):
        '''
        Function to mesh a part using reduced order C3D8R elements with
        enhanced hourglass control and distortion control.

        Args:
            part (Abaqus Part Instance): the part to be meshed.
            mesh_size (float): size of elements in the FE mesh.

        Returns:
            None
        '''
        elem_library = EXPLICIT
        region = part.Set(
            cells=part.cells, name='{}-Cells'.format(part.name))
        elemType1 = mesh.ElemType(
            elemCode=C3D8R, elemLibrary=elem_library, lengthRatio=0.10,
            kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
            hourglassControl=ENHANCED, distortionControl=ON)
        elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=elem_library)
        elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=elem_library)
        part.setElementType(
            regions=region,
            elemTypes=(elemType1, elemType2, elemType3))
        part.setMeshControls(
            regions=part.cells, algorithm=ADVANCING_FRONT,
            elemShape=HEX_DOMINATED, technique=SWEEP)
        part.seedPart(
            size=self.solid_mesh, deviationFactor=0.1, minSizeFactor=0.1)
        part.generateMesh()


if __name__ == '__main__':
    # Material parameters
    tape_angles = (0, 45, 90, -45)   # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.205  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.0683333  # ratio of undulation amplitude to length
    number_of_plies = 8  # symmetric 4 ply laminate

    # RVE dimensions
    x_min = 0.0
    x_max = 25.0
    y_min = 0.0
    y_max = 25.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = AP_PLY_3D('TestModel')
    mdl.symmetry_mode = 'BC'
    mdl.solid_mesh = 1.5  # set mesh size for testing purposes
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies, x_shift=-12.5, y_shift=-12.5)
    mdl.create_materials()
    paths_f, angles_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_parts(grid_f, mesh=True)
    mdl.create_sequences()
    mdl.create_cohesive_interactions(tie=True)
    mdl.create_job()
    mdl.save_model()
