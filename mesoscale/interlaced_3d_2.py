'''
This module can be used to create 2D (shell) models of interlaced laminates.
The Interlaced2D class inherits from the InterlacedModel class. The child
class contains methods to create parts, partition them, and assign material
properties to the partitioned regions.

(c) Rutger Kok, 21/06/2021
'''
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\Github\\interlaced_model_creation\\mesoscale')
from abaqus import *
from abaqusConstants import *
from caeModules import *
import mesh
import regionToolset
from shapely.geometry import Polygon
from shapely.ops import unary_union
from interlaced_model import InterlacedModel
import numpy as np


class Interlaced3D(InterlacedModel):
    def __init__(self, model_name, dir_name=None):
        InterlacedModel.__init__(self, model_name, dir_name)

    def create_parts(self, grid):
        '''
        This function is used to create parts using geometric info from
        Shapely.
        Args:
            grid (dict): dictionary containing Shapely polygons

        Returns:
            None
        '''
        region_types = ['Tape', 'Resin', 'Undulation']
        for n in range(1, len(self.t_angles) + 1):
            instances = []  # instances in layer
            mirror_instances = []  # instances in mirrored layer
            name = 'Layer-{}'.format(n)  # layer name
            mirror_name = name + '-Mirror'
            # define dictionary to store material orientations
            self.mat_orient = {name: {}, mirror_name: {}}
            for obj_type in region_types:
                for angle in self.t_angles:
                    # determine all objects of same type in a single layer
                    same_type = [
                        obj for obj in grid.values()
                        if obj.object_type[n - 1] == obj_type
                        and obj.angle[n - 1] == angle]
                    merged_obj = unary_union(same_type)
                    # split merged_obj into  its constituent Polygons
                    polygons = self.split_multipolygon(merged_obj)
                    # iterate over polygons to create parts
                    for i, poly in enumerate(polygons):
                        parts = []
                        try:
                            part_id = '{}-{}-{}-{}'.format(
                                obj_type, angle, n, i)
                            part = self.create_3d_part(
                                poly, part_id, obj_type, n, angle)
                            parts.append((
                                part_id, part,
                                poly.representative_point().coords[0]))
                        except:
                            individual_polys = [
                                obj for obj in grid.values()
                                if obj.within(
                                    poly.buffer(1.0e-6, join_style=2))]
                            for j, sub_poly in enumerate(individual_polys):
                                part_id = '{}-{}-{}-{}-{}'.format(
                                    obj_type, angle, n, i, j)
                                part = self.create_3d_part(
                                    sub_poly, part_id, obj_type, n,
                                    angle)
                                parts.append((
                                    part_id, part,
                                    sub_poly.representative_point().coords[0]))
                        # create part instances
                        for part_id, part, centroid in parts:
                            instance, mirror_instance = self.create_instances(
                                part, part_id, angle, centroid, n)
                            instances.append(instance)
                            if mirror_instance:
                                mirror_instances.append(mirror_instance)
            # merge part instances in a single layer into one part
            self.merge_layer(name, instances)
            if mirror_instances:
                self.merge_layer(
                    mirror_name, mirror_instances)

    def merge_layer(self, layer_name, instances):
        self.assembly.InstanceFromBooleanMerge(
            name=layer_name, keepIntersections=ON, instances=instances,
            domain=GEOMETRY, originalInstances=DELETE)
        merged_part = self.model.parts[layer_name]
        self.assign_orientations(merged_part, layer_name)
        self.mesh_part_3d(merged_part, self.mesh_size)

    def assign_orientations(self, merged_part, layer_name):
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
        sequence_thickness = len(self.t_angles) * self.t_thickness
        for sequence in range(self.sequences - 1):
            instance_names = [
                instance for instance in self.instances.keys()
                if 'Mirror' not in instance]
            self.assembly.translate(
                instanceList=instance_names,
                vector=(0.0, 0.0, sequence * sequence_thickness))
            if self.symmetry_mode == 'GEOM':
                mirror_instance_names = [
                    instance for instance in self.assembly.instances.keys()
                    if 'Mirror' in instance]
                self.assembly.translate(
                    instanceList=mirror_instance_names,
                    vector=(0.0, 0.0, -1.0 * sequence * sequence_thickness))

    def create_cohesive_interactions(self, tie=False):
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
            if self.damage_mode == 'ON':
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

    def create_3d_part(self, obj, obj_id, obj_type, obj_layer, obj_angle):
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
            part.BaseSolidExtrude(sketch=sketch, depth=self.t_thickness)
        cell_region = regionToolset.Region(cells=part.cells[:])
        if self.damage_mode == 'ON':
            mat_type = 'Damage'
        else:
            mat_type = 'Elastic'
        if obj_type == 'Tape':
            section_name = 'Tape-{}'.format(mat_type)
        elif obj_type == 'Resin':
            section_name = 'Resin-Tape-Elastic'  # always elastic
        elif obj_type == 'Undulation':
            section_name = 'Undulation-{}-(0, Resin)'.format(mat_type)
        part.SectionAssignment(
            region=cell_region, sectionName=section_name,
            offsetType=MIDDLE_SURFACE, thicknessAssignment=FROM_SECTION,
            offset=0.0, offsetField='')
        return part

    def create_instances(self, part, obj_id, obj_angle, rep_point, obj_layer):
        # create the part instance
        z_offset = (obj_layer - 1) * self.t_thickness
        instance = self.assembly.Instance(
            name=obj_id, part=part, dependent=ON)
        self.assembly.translate(
            instanceList=(obj_id, ), vector=(0.0, 0.0, z_offset))
        z_coord = self.t_thickness * (2 * obj_layer - 1) / 2.0
        self.mat_orient['Layer-{}'.format(obj_layer)].setdefault(
            obj_angle, []).append((rep_point[0], rep_point[1], z_coord))
        # if symmetry mode is set to geometric rather than boundary condition
        # then mirrored part instances are created
        if self.symmetry_mode == 'GEOM':
            # create the symmetric part instance
            mirror_instance = self.assembly.Instance(
                name='Mirror-' + obj_id, part=part, dependent=ON)
            translate_vector = (0.0, 0.0, -z_offset - self.t_thickness)
            self.assembly.translate(
                instanceList=('Mirror-' + obj_id, ),
                vector=translate_vector)
            self.mat_orient['Layer-{}-Mirror'.format(obj_layer)].setdefault(
                obj_angle, []).append((rep_point[0], rep_point[1], -z_coord))
            return [instance, mirror_instance]
        return [instance, None]

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
        # first try simplify using convex hull
        ch_polygon = poly_in.convex_hull
        if ch_polygon.area == poly_in.area:
            return ch_polygon
        else:
            # if convex hull oversimplifies polygon then adjust vertices by
            # their angles
            ext_poly_coords = poly_in.exterior.coords[:]
            vector_rep = np.diff(ext_poly_coords, axis=0)
            angles_list = []
            for i in range(0, len(vector_rep) - 1):
                angles_list.append(
                    np.abs(self.get_angles(vector_rep[i], vector_rep[i + 1])))
            # get mask satisfying tolerance
            thresh_vals_by_deg = np.where(np.array(angles_list) > deg_tol)
            # gotta be a better way to do this next part
            # sandwich betweens first and last points
            new_idx = [0] + (thresh_vals_by_deg[0] + 1).tolist() + [0]
            new_vertices = [ext_poly_coords[idx] for idx in new_idx]
            clean_polygon = Polygon(new_vertices)
            return clean_polygon

    def mesh_part_3d(self, part, mesh_size):
        '''
        Function to mesh a part using reduced order C3D8R elements with
        enhanced hourglass control and distortion control.

        Args:
            part (Abaqus Part Instance): the part to be meshed.
            mesh_size (float): size of elements in the FE mesh.

        Returns:
            None
        '''
        if self.damage_mode == 'ON':
            elem_library = STANDARD
        else:
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
            regions=region, elemTypes=(elemType1, elemType2, elemType3))
        min_volume = (self.u_width ** 2.0) * self.t_thickness
        i, j = 0, 0
        for cell in part.cells:
            # calculate cell volume and shape to assign wedge element type to
            # small volume regions. Remaining regions are meshed using
            # hex elements
            cell_volume = cell.getSize(printResults=False)
            if cell_volume <= min_volume:
                if i == 0:
                    wedge_cells = part.cells.findAt(cell.pointOn, )
                    i = 1
                else:
                    wedge_cells += part.cells.findAt(cell.pointOn, )
            else:
                if j == 0:
                    hex_cells = part.cells.findAt(cell.pointOn, )
                    j = 1
                else:
                    hex_cells += part.cells.findAt(cell.pointOn, )
        part.setMeshControls(
            regions=hex_cells, algorithm=MEDIAL_AXIS,
            elemShape=HEX_DOMINATED, technique=SWEEP)
        if i > 0:
            part.setMeshControls(
                regions=wedge_cells, elemShape=WEDGE)    
        part.seedPart(
            size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
        part.generateMesh()


if __name__ == '__main__':
    # Material parameters
    tape_angles = (0, 90)   # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 4 ply laminate

    # Simulation parameters
    damage_mode = 'OFF'

    # RVE dimensions
    x_min = -20.0
    x_max = 20.0
    y_min = 20.0
    y_max = 20.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = Interlaced3D('TestModel')
    mdl.damage_mode = damage_mode
    mdl.symmetry_mode = 'BC'
    mdl.mesh_size = 0.5  # set mesh size for testing purposes
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    paths_f, angles_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_parts(grid_f)
    mdl.create_cohesive_interactions()
    mdl.create_job()
    mdl.save_model()

