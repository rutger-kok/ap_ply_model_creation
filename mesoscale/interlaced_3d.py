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
from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from math import cos, radians, tan
from interlaced_model import InterlacedModel


class Interlaced3D(InterlacedModel):
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
                        faces=part.faces)
                else:
                    part.PartitionCellByDatumPlane(
                        datumPlane=part.datums[datum_plane],
                        cells=part.cells)
            except:
                # bare except to catch feature creation error which occurs
                # when the partitioning plane does not intercept the part
                 continue

    def partition_through_thickness(self):
        for n in range(1, self.l_plies):
            datum_plane_id = self.specimen_part.DatumPlaneByPrincipalPlane(
                principalPlane=XYPLANE, offset=self.t_thickness * n).id
            self.specimen_part.PartitionCellByDatumPlane(
                datumPlane=self.specimen_part.datums[datum_plane_id],
                cells=self.specimen_part.cells)

    def create_part(self, x_min, y_min, x_max, y_max, part_grid):
        '''
        Function to create a 3d rectangular specimen part.

        Args:
            x_min (float): minimum x-coordinate of the specimen.
            x_max (float): maximum x-coordinate of the specimen.
            y_min (float): minimum y-coordinate of the specimen.
            y_max (float): maximum y-coordinate of the specimen.
            part_grid (dict): grid dictionary that defines the properties of
                the interlaced regions.

        Returns:
            None
        '''
        part_sketch = self.model.ConstrainedSketch(
            name='Part Sketch', sheetSize=200.0)
        part_sketch.rectangle(point1=(x_min, y_min), point2=(x_max, y_max))
        self.specimen_part = self.model.Part(
            name='Specimen', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.specimen_part.BaseSolidExtrude(
            sketch=part_sketch, depth=self.l_thickness)
        self.specimen_faces = self.specimen_part.faces
        self.assembly.Instance(
            name='Specimen Instance', part=self.specimen_part, dependent=ON)
        # partition the part into subregions
        self.partition_part(self.specimen_part, face=False)
        # parition the part into plies through the thickness
        self.partition_through_thickness()
        # assign properties to each of the regions
        self.assign_properties(part_grid)
        # mesh the part
        # self.mesh_part_3d(self.specimen_part, self.mesh_size)

    def assign_properties(self, grid):
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
        cells = self.specimen_part.cells
        if self.damage_mode == 'ON':
            mat_type = 'Damage'
        else:
            mat_type = 'Elastic'
        region_dict = {}
        orient_dict = {}
        for obj_id, obj in grid.iteritems():
            obj_centroid = obj.centroid.coords
            for n in range(1, len(self.t_angles) + 1):
                region_coords = (
                    obj_centroid[0][0], obj_centroid[0][1],
                    self.t_thickness * (2 * n - 1) / 2.0)
                obj_type = obj.object_type[n - 1]
                # check object type and assign properties accordingly
                if obj_type == 'Tape':
                    section_name = 'Tape-{}'.format(mat_type)
                    mat_orient = obj.angle[n - 1]
                elif obj_type == 'Resin':
                    mat_orient = obj.angle[n - 1]
                    section_name = 'Resin-Tape-Elastic'  # always elastic
                elif obj_type == 'Undulation':
                    mat_orient = obj.angle[n - 1]
                    section_name = 'Undulation-{}-(0, Resin)'.format(mat_type)
                # if region key does not exist in dictionary create it, else
                # append the centroid of the new region to the dictionary entry
                region_dict.setdefault(section_name, []).append(region_coords)
                orient_dict.setdefault(mat_orient, []).append(region_coords)

        for name, regions in region_dict.iteritems():
            for i, centroid in enumerate(regions):
                if i == 0:
                    selected_cells = cells.findAt((centroid, ))
                else:
                    selected_cells += cells.findAt((centroid, ))
            # assign section to region
            region = self.specimen_part.Set(
                cells=selected_cells, name='{}-Cells'.format(name))
            self.specimen_part.SectionAssignment(
                region=region, sectionName=name, offsetType=MIDDLE_SURFACE,
                offset=0.0, offsetField='', thicknessAssignment=FROM_SECTION)

        for orientation, regions in orient_dict.iteritems():
            for i, centroid in enumerate(regions):
                if i == 0:
                    selected_cells = cells.findAt((centroid, ))
                else:
                    selected_cells += cells.findAt((centroid, ))
            # assign orientation to region
            region = self.specimen_part.Set(
                cells=selected_cells, name='Cells-{}'.format(int(orientation)))
            self.specimen_part.MaterialOrientation(
                region=region, orientationType=SYSTEM, localCsys=None,
                stackDirection=STACK_3, additionalRotationField='',
                angle=orientation, additionalRotationType=ROTATION_ANGLE,
                axis=AXIS_3, fieldName='')

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
            element_library = STANDARD
        else:
            element_library = EXPLICIT
        elem_type_1 = mesh.ElemType(
            elemCode=C3D8R, elemLibrary=element_library, lengthRatio=0.1,
            kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF,
            hourglassControl=ENHANCED, distortionControl=ON)
        elem_type_2 = mesh.ElemType(elemCode=C3D6, elemLibrary=element_library)
        elem_type_3 = mesh.ElemType(elemCode=C3D4, elemLibrary=element_library)
        part_cells = part.Set(
            cells=part.cells[:], name='{} Cells'.format(part.name))
        part.setElementType(
            regions=part_cells,
            elemTypes=(elem_type_1, elem_type_2, elem_type_3))
        part.seedPart(
            size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
        part.generateMesh()


if __name__ == '__main__':
    # Material parameters
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 2  # symmetric 8 ply laminate

    # Simulation parameters
    damage_mode = 'OFF'

    # RVE dimensions
    x_min = -50.0
    x_max = 50.0
    y_min = -75.0
    y_max = 75.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = Interlaced3D('TestModel')
    mdl.damage_mode = damage_mode
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    mdl.mesh_size = 0.5  # set mesh size for testing purposes
    paths_f, angles_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_part(x_min, y_min, x_max, y_max, grid_f)
    mdl.create_job()
    mdl.save_model()
