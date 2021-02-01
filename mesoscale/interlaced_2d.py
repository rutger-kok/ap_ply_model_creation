'''
This module can be used to create 2D (shell) models of interlaced laminates.
The Interlaced2D class inherits from the InterlacedModel class. The child
class contains methods to create parts, partition them, and assign material
properties to the partitioned regions.

(c) Rutger Kok, 25/11/2020
'''
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\Github\\interlaced_model_creation\\mesoscale')
from abaqus import *
from abaqusConstants import *
from caeModules import *
import mesh
import regionToolset
from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from math import cos, radians, tan
from interlaced_model import InterlacedModel


class Interlaced2D(InterlacedModel):
    def __init__(self, model_name, dir_name=None):
        InterlacedModel.__init__(self, model_name, dir_name)

    def partition_part(self, part, face=True):
        '''
        Function to partition the part into a grid of polygons depending on the
        tape angles of the laminate.

        Args:
            part (Abaqus Part Instance): part to be partitioned.
            face (boolean): True if partitioning faces (2D parts), False if
                partitioning cells.

        Returns:
            None
        '''
        # define mirror point (used to mirror the Polygon boundaries)
        mirror_point = Point([(0.0, 0.0), (0.0, 0.0)])

        x, y = zip(*self.specimen_size.exterior.coords)
        x_max = max(x)  # determine max x-value for offsets
        y_max = max(y)

        partition_lines = []
        uw = self.u_width
        for a, w in zip(self.t_angles, self.t_widths):
            if a == 90:
                offset = w
                max_offset = x_max - w / 2.0
                number_offsets = int(max_offset / offset)
                offset_list = [offset * k for k in range(number_offsets + 1)]
                for ofs in offset_list:
                    line_1_coords = [((w / 2.0) - uw + ofs, 100.0),
                                     ((w / 2.0) - uw + ofs, -100.0)]
                    line_2_coords = [((w / 2.0) + ofs, 100.0),
                                     ((w / 2.0) + ofs, -100.0)]
                    line_3_coords = [((w / 2.0) + ofs + uw, 100.0),
                                     ((w / 2.0) + ofs + uw, -100.0)]
                    line_1 = LineString(line_1_coords)
                    line_2 = LineString(line_2_coords)
                    line_3 = LineString(line_3_coords)
                    partition_lines.extend([line_1, line_2, line_3])
                reflected_lines = [scale(line, xfact=-1, origin=mirror_point)
                                   for line in partition_lines]
                partition_lines = partition_lines + reflected_lines

            else:
                offset = w / cos(radians(a))
                max_offset = (y_max - (w / 2.0) * cos(radians(a))
                              + x_max * tan(radians(a)))
                number_offsets = int(max_offset / offset)
                offset_list = [offset * k for k in range(number_offsets + 20)]
                for ofs in offset_list:
                    line_1_coords = [(-100.0, (w / 2.0) - uw + ofs),
                                     (100.0, (w / 2.0) - uw + ofs)]
                    line_2_coords = [(-100.0, (w / 2.0) + ofs),
                                     (100.0, (w / 2.0) + ofs)]
                    line_3_coords = [(-100.0, (w / 2.0) + ofs + uw),
                                     (100.0, (w / 2.0) + ofs + uw)]
                    rotation_point = Point([(0.0, ofs), (0.0, ofs)])
                    line_1 = rotate(
                        LineString(line_1_coords), a, rotation_point)
                    line_2 = rotate(
                        LineString(line_2_coords), a, rotation_point)
                    line_3 = rotate(
                        LineString(line_3_coords), a, rotation_point)
                    partition_lines.extend([line_1, line_2, line_3])
                reflected_lines = [rotate(line, 180.0, origin=mirror_point)
                                   for line in partition_lines]
                partition_lines = partition_lines + reflected_lines

        for line in partition_lines:
            line_coords = [list(tup) + [0.0] for tup in list(line.coords)]
            last_coord = [[line_coords[-1][0], line_coords[-1][1], 1.0]]
            line_coords_3d = line_coords + last_coord
            datum_point_ids = []
            for coordinates in line_coords_3d:
                datum_point_ids.append(
                    part.DatumPointByCoordinate(
                        coords=coordinates).id)
            datum_point_1 = part.datums[datum_point_ids[0]]
            datum_point_2 = part.datums[datum_point_ids[1]]
            datum_point_3 = part.datums[datum_point_ids[2]]
            datum_plane = part.DatumPlaneByThreePoints(
                point1=datum_point_1, point2=datum_point_2,
                point3=datum_point_3).id
            try:
                if face:
                    part.PartitionFaceByDatumPlane(
                        datumPlane=part.datums[datum_plane],
                        faces=self.specimen_faces)
                else:
                    part.PartitionCellByDatumPlane(
                        datumPlane=part.datums[datum_plane],
                        cells=part.cells)
            except:
                continue

    def create_part(self, x_min, y_min, x_max, y_max, part_grid):
        '''
        Function to create a 2d (shell) rectangular specimen part.

        Args:
            x_min (float): minimum x-coordinate of the specimen.
            x_max (float): maximum x-coordinate of the specimen.
            y_min (float): minimum y-coordinate of the specimen.
            y_max (float): maximum y-coordinate of the specimen.
            part_grid (dict): dictionary of grids that define properties of
                interlaced regions.

        Returns:
            None
        '''
        part_sketch = self.model.ConstrainedSketch(
            name='Part Sketch', sheetSize=200.0)
        part_sketch.rectangle(point1=(x_min, y_min), point2=(x_max, y_max))
        self.specimen_part = self.model.Part(
            name='Specimen', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.specimen_part.BaseShell(sketch=part_sketch)
        self.specimen_faces = self.specimen_part.faces
        self.assembly.Instance(
            name='Specimen Instance', part=self.specimen_part, dependent=ON)

        self.partition_part(self.specimen_part)
        self.mesh_part_2d(self.specimen_part, self.mesh_size)
        self.assign_properties(self.specimen_part, part_grid)

    def assign_properties(self, part, part_grid):
        '''
        Function to assign material properties to the partitioned specimen.

        Args:
            part (Abaqus Part Instance): part for property assignment.
            part_grid (dict): Nested dictionary defining the grid of
                tape/undulation objects.

        Returns:
            None

        TODO:
            Change property assignment to assign to all regions of the same
                type at once (reducing the number of composite layups)
        '''
        faces = part.faces
        for obj_id, obj in part_grid[1].iteritems():
            obj_centroid = obj.centroid.coords
            selected_face = faces.findAt(
                ((obj_centroid[0][0], obj_centroid[0][1], 0.0), ))
            region = regionToolset.Region(faces=selected_face)
            composite_layup = part.CompositeLayup(
                name='CompositeLayup-{}'.format(obj_id), description='',
                elementType=SHELL, offsetType=MIDDLE_SURFACE, symmetric=False,
                thicknessAssignment=FROM_SECTION)
            composite_layup.Section(
                preIntegrate=OFF, integrationRule=SIMPSON, useDensity=OFF,
                thicknessType=UNIFORM, poissonDefinition=DEFAULT,
                temperature=GRADIENT)
            composite_layup.ReferenceOrientation(
                orientationType=GLOBAL, localCsys=None, axis=AXIS_3,
                fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0)
            for layer_number in range(1, len(part_grid) + 1):
                layer_obj = part_grid[layer_number][obj_id]
                obj_type = layer_obj.object_type
                obj_angle = layer_obj.angle[-1]
                if obj_type == 'Tape':
                    material_name = 'Tape-Elastic'
                elif obj_type == 'Undulation':
                    if len(layer_obj.angle) == 1:
                        material_name = 'Undulation-Elastic-Resin'
                    else:
                        interface_angle = abs(
                            layer_obj.angle[0] - layer_obj.angle[1])
                        section_angle = (0, interface_angle)
                        material_name = 'Undulation-Elastic-{}'.format(
                            section_angle)
                else:
                    material_name = 'Resin-Elastic'
                name_of_ply = 'Ply-{}'.format(layer_number)
                composite_layup.CompositePly(
                    suppressed=False, plyName=name_of_ply, region=region,
                    material=material_name, thicknessType=SPECIFY_THICKNESS,
                    thickness=0.2, orientationType=SPECIFY_ORIENT, angle=0.0,
                    orientationValue=obj_angle, additionalRotationField='',
                    additionalRotationType=ROTATION_NONE, numIntPoints=3,
                    axis=AXIS_3)

        # # find all combinations of angles in a layer
        # for layer_id, layer_grid in part_grid.iteritems():
        #     angle_sets = set()
        #     for obj_id, obj in layer_grid.iteritems():
        #         angle_sets.add(tuple(obj.angle))

        #     # find faces with the same angle
        #     for i, angle in enumerate(angle_sets):
        #         obj_by_angle = [obj for obj in layer_grid.values()
        #                         if tuple(obj.angle) == angle]
        #         # now split regions with same angle into lists by object type
        #         object_types = ('Tape', 'Resin', 'Undulation')
        #         obj_by_type = {
        #             region_type: [obj.centroid.coords for obj in obj_by_angle
        #                           if obj.object_type == region_type]
        #             for region_type in object_types}
        #         # iterate over sorted objects to assign properties
        #         for region_type, centroids in obj_by_type.iteritems():
        #             if region_type == 'Tape':
        #                 section_name = 'Tape-Elastic'
        #             elif region_type == 'Undulation':
        #                 if len(angle) == 1:
        #                     section_name = 'Undulation-Elastic-Resin'
        #                 else:
        #                     interface_angle = abs(angle[0] - angle[1])
        #                     section_angle = (0, interface_angle)
        #                     section_name = 'Undulation-Elastic-{}'.format(
        #                         section_angle)
        #             else:
        #                 section_name = 'Resin-Elastic'
        #             if centroids:
        #                 selected_faces = self.specimen_faces.findAt(
        #                     ((centroids[0][0][0], centroids[0][0][1], 0.0), ))
        #                 for n in range(1, len(centroids)):
        #                     selected_faces += self.specimen_faces.findAt(
        #                         ((centroids[n][0][0], centroids[n][0][1], 0.0), ))
        #                 region = regionToolset.Region(faces=selected_faces)
        #                 composite_layup = self.specimen_part.CompositeLayup(
        #                     name='CompositeLayup-{}'.format(i), description='',
        #                     elementType=SHELL, offsetType=MIDDLE_SURFACE,
        #                     symmetric=False, thicknessAssignment=FROM_SECTION)
        #                 composite_layup.Section(
        #                     preIntegrate=OFF, integrationRule=SIMPSON,
        #                     thicknessType=UNIFORM, poissonDefinition=DEFAULT,
        #                     temperature=GRADIENT, useDensity=OFF)
        #                 composite_layup.ReferenceOrientation(
        #                     orientationType=GLOBAL, fieldName='', angle=0.0,
        #                     axis=AXIS_3, additionalRotationType=ROTATION_NONE,
        #                     localCsys=None)
        #                 for j in range(self.sequences):
        #                     name_of_ply = 'Ply-{}-{}'.format(j, layer_id)
        #                     composite_layup.CompositePly(
        #                         region=region, material=section_name,
        #                         thicknessType=SPECIFY_THICKNESS, numIntPoints=3,
        #                         thickness=self.t_thickness, suppressed=False,
        #                         orientationType=SPECIFY_ORIENT, angle=0.0,
        #                         orientationValue=0.0, plyName=name_of_ply,
        #                         additionalRotationType=ROTATION_NONE, axis=AXIS_3,
        #                         additionalRotationField='')

    def mesh_part_2d(self, part, mesh_size):
        '''
        Function used to mesh the 2D part using shell elements (with enhanced
        hourglass control)

        Args:
            part (Abaqus Part Instance): part to assign mesh to.
            mesh_size (float): size of elements in FE mesh.

        Returns:
            None
        '''
        face_region = regionToolset.Region(faces=part.faces[:])
        elemType1 = mesh.ElemType(
            elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF,
            hourglassControl=ENHANCED)
        elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
        part.setElementType(
            regions=face_region, elemTypes=(elemType1, elemType2))
        part.seedPart(size=mesh_size)
        part.generateMesh()


if __name__ == '__main__':
    # Material parameters
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 15.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 8 ply laminate

    # RVE dimensions
    x_min = -50.0
    x_max = 50.0
    y_min = -75.0
    y_max = 75.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = Interlaced2D('TestModel')
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    mdl.mesh_size = 0.5  # set mesh size for testing purposes
    paths_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_part(x_min, y_min, x_max, y_max, grid_f)
    mdl.create_job()
    mdl.save_model()
