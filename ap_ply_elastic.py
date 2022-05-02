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
sys.path.append('C:\\Python27\\Lib\\site-packages')
sys.path.append('C:\\Github\\ap_ply_model_creation')
from abaqus import *
from abaqusConstants import *
from caeModules import *
import mesh
import regionToolset
from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from shapely.ops import unary_union, polygonize
from math import cos, radians, tan, sin
from ap_ply_model import AP_PLY_Model


class AP_PLY_Elastic(AP_PLY_Model):
    def __init__(self, model_name, dir_name=None):
        AP_PLY_Model.__init__(self, model_name, dir_name)

    def create_grid_elastic(self, region):
        '''
        Method to create a grid of Shapely polygons depending on
        the tape angles of the laminate (note: no undulation regions!).
        The region region (XXX mm x XXX mm) is separated into Shapely Polygons.
        '''
        region_buffered = region.buffer(
            2.0 * self.u_width, join_style=2)
        x, y = zip(*region_buffered.exterior.coords)
        # NOTE: these max values are different from the x_min etc. values
        x_max = max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)

        x_len = max(abs(self.x_0 - x_max), abs(self.x_0 - x_min))
        y_len = max(abs(self.y_0 - y_max), abs(self.y_0 - y_min))

        tapes = zip(self.t_angles, self.t_widths)
        partition_lines = []
        for (a, w) in tapes:
            if a == 90 or a == -90.0:
                n_offsets = int(x_len / w) + 1
                offset_list = [w * k for k in range(-n_offsets, n_offsets)]
                for ofs in offset_list:
                    x_2 = self.x_0 + (w / 2.0) + ofs
                    coords = [(x_2, y_max), (x_2, y_min)]
                    line = LineString(coords)
                    try:
                        line_trim = line.intersection(region_buffered)
                        if line_trim:
                            partition_lines.append(line_trim)
                    except NotImplementedError:
                        continue
            else:
                a_rad = radians(a)
                offset = w / cos(a_rad)
                max_offset = (y_len + abs(x_len * tan(a_rad)))
                n_offsets = int(max_offset / offset) + 1
                offset_list = [
                    offset * k for k in range(-n_offsets, n_offsets)]
                for ofs in offset_list:
                    y_2 = self.y_0 + (w / 2.0) + ofs
                    line = self.create_line(
                        y_2, a_rad, ofs, x_min, x_max, region_buffered)
                    if line:
                        partition_lines.append(line)

        # collection of individual linestrings for splitting in a list and add
        # the polygon lines to it.
        partition_lines.append(region.boundary)
        partition_lines.append(region_buffered.boundary)
        border_lines = unary_union(partition_lines)
        decomposition = polygonize(border_lines)
        trimmed = [polygon for polygon in decomposition
                   if polygon.within(region)]
        polygon_dict = {idx: polygon for (idx, polygon) in enumerate(trimmed)}
        # remove sliver polygons
        min_area = 1.0e-4
        clean_dict = {
            idx: polygon for (idx, polygon) in polygon_dict.iteritems()
            if polygon.area > min_area}
        return clean_dict, partition_lines

    def partition_part(self, part, partition_lines, cells=True):
        '''
        Method to parition an Abaqus (2D shell or 3D solid) part along the
        lines specified in the partition_lines list.
        '''
        # Partition part
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
                if cells:
                    part.PartitionCellByDatumPlane(
                        datumPlane=part.datums[datum_plane],
                        cells=part.cells)
                else:
                    part.PartitionFaceByDatumPlane(
                        datumPlane=part.datums[datum_plane],
                        faces=part.faces)
            except:
                continue

    def partition_through_thickness(self, part):
        '''
        Method to partition a solid part through the thickness into separate
        plies.
        '''
        # Partition part through the thickness
        if self.symmetry_mode == 'GEOM':
            offsets = [
                0.0 + i * self.t_thickness for i
                in range(1, self.l_plies)]
        else:
            offsets = [
                0.0 + i * self.t_thickness for i
                in range(1, self.l_plies / 2)]
        for offset in offsets:
            dplane = part.DatumPlaneByPrincipalPlane(
                    principalPlane=XYPLANE, offset=offset).id
            try:
                part.PartitionCellByDatumPlane(
                    datumPlane=part.datums[dplane], cells=part.cells)
            except:
                continue

    def create_part_3d(self, region):
        '''
        Method to create a 3D solid part which can be partitioned into
        regions according to the AP-PLY architecture.
        '''
        x, y = zip(*region.exterior.coords)
        x_max = max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)

        x_mid = (x_max + x_min) / 2.0
        y_mid = (y_max + y_min) / 2.0
        # define origin point (where pattern starts)
        self.x_0 = x_mid + self.x_shift
        self.y_0 = y_mid + self.y_shift

        if self.symmetry_mode == 'BC':
            z_min = 0.0
            z_max = self.l_thickness / 2.0
        else:
            z_min = -self.l_thickness / 2.0
            z_max = self.l_thickness / 2.0
        sketch = self.model.ConstrainedSketch(
            name='Sketch', sheetSize=200.0)
        sketch.rectangle(
            point1=(x_min, y_min), point2=(x_max, y_max))
        # Create a 3D deformable part
        part = self.model.Part(
            name='Specimen', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseSolidExtrude(sketch=sketch, depth=(z_max - z_min))
        self.assembly.Instance(
            name='Specimen Instance', part=part, dependent=ON)

        return part

    def make_pass_elastic(self, grid, region, angle=0, pass_coords=None):
        '''
        Method used to virtually place an AFP machine pass.
        NOTE: a pass should be specified using its coordinates when horizontal
        and then rotated. Do not rotate outside of the function
        '''
        # default tape dimensions (used unless other coords are supplied)
        x, y = zip(*region.exterior.coords)
        x_max = max(x)
        x_min = min(x)
        if pass_coords is None:
            w = self.t_widths[0]
            pass_coords = ([
                (x_min, self.y_0 + (w / 2)),
                (x_max, self.y_0 + (w / 2)),
                (x_max, self.y_0 + (-w / 2)),
                (x_min, self.y_0 + (-w / 2))])

        # define tape boundaries
        tol = 0.001  # tolerance for buffering
        pass_bounds = scale(Polygon(pass_coords), xfact=5.0)
        pass_rotated_bounds = rotate(
            pass_bounds, angle).buffer(tol, join_style=2)
        # identify which Polygons from the create_grid_2d function are within
        # the boundaries of the rotated tape
        if pass_rotated_bounds.intersects(region):
            pass_path = [(obj_id, obj) for (obj_id, obj)
                         in grid.iteritems()
                         if obj.within(pass_rotated_bounds)]
            for i in range(len(pass_path)):
                objID, obj = pass_path[i]
                if not hasattr(obj, 'angle'):
                    setattr(obj, 'angle', [angle, ])
                else:
                    obj.angle.append(angle)

    def place_tows_elastic(self, grid, region):
        '''
        Method which defines the coordinates of the tow paths to be placed
        based on the AP-PLY laminate parameters.
        '''
        x, y = zip(*region.exterior.coords)
        x_max = max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)
        x_len = max(abs(self.x_0 - x_max), abs(self.x_0 - x_min))
        y_len = max(abs(self.y_0 - y_max), abs(self.y_0 - y_min))

        tape_props = zip(self.t_angles, self.t_widths)
        s = self.t_spacing  # reassign t_spacing for brevity
        for n_shift in range(s + 1):
            for (a, w) in tape_props:
                a_rad = radians(a)
                if a == 90.0 or a == -90.0:
                    offset = w
                    spaced_offset = (1 + s) * w
                    # largest offset within region
                    max_offset = x_len
                else:
                    offset = w / cos(a_rad)
                    spaced_offset = ((1 + s) * w) / cos(a_rad)
                    max_offset = (y_len + (w / 2.0) + abs(x_len * tan(a_rad)))
                n_offsets = int(max_offset / spaced_offset) + 5
                offset_list = [
                    spaced_offset * k for k in range(-n_offsets, n_offsets)]
                shift = n_shift * offset
                for ofs in offset_list:
                    if a == 90.0 or a == -90.0:
                        x_1 = self.x_0 - x_len + ofs + shift
                        x_2 = self.x_0 + x_len + ofs + shift
                        y_1 = w / 2.0
                        y_2 = -w / 2.0
                    else:
                        x_1 = self.x_0 - x_len
                        x_2 = self.x_0 + x_len
                        y_1 = self.y_0 + (w / 2.0) + ofs + shift
                        y_2 = self.y_0 + (-w / 2.0) + ofs + shift
                    tape_coords = ([(x_1, y_1), (x_2, y_1),
                                    (x_2, y_2), (x_1, y_2)])
                    self.make_pass_elastic(
                        grid, region, angle=a, pass_coords=tape_coords)

    def create_line(self, y_offset, a_rad, offset, x_min, x_max, region):
        '''
        Method used to create the lines required to partition the laminate
        geometry.
        '''
        y_delta = sin(a_rad) * (y_offset - (self.y_0 + offset))
        x_neg = ((((x_min - self.x_0) + y_delta) / cos(a_rad)) + self.x_0)
        x_pos = ((((x_max - self.x_0) + y_delta) / cos(a_rad)) + self.x_0)
        coords = [(x_neg, y_offset), (x_pos, y_offset)]
        line = rotate(
            LineString(coords), a_rad, Point([self.x_0, self.y_0 + offset]),
            use_radians=True)
        scaled_line = scale(line, xfact=2.0, yfact=2.0)
        try:
            line_trim = scaled_line.intersection(region)
            return line_trim
        except NotImplementedError:
            return None

    def create_part_2d(self, shell_region):
        '''
        Method to create a 2d (shell) rectangular shell_region part.
        '''
        self.shell_region = shell_region
        x, y = zip(*self.shell_region.exterior.coords)
        # NOTE: these max values are different from the x_min etc. values
        x_max = max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)

        x_mid = (x_max + x_min) / 2.0
        y_mid = (y_max + y_min) / 2.0
        # define origin point (where pattern starts)
        self.x_0 = x_mid + self.x_shift
        self.y_0 = y_mid + self.y_shift

        part_sketch = self.model.ConstrainedSketch(
            name='Part Sketch', sheetSize=200.0)
        part_sketch.rectangle(
            point1=(x_min, y_min), point2=(x_max, y_max))
        part = self.model.Part(
            name='Specimen', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseShell(sketch=part_sketch)
        self.assembly.Instance(
            name='Specimen Instance', part=part, dependent=ON)
        return part

    def group_by_angle_3d(self, grid):
        '''
        Method used to group ojects of the same angle to allow for more
        efficient material orientation assignment.
        '''
        if self.symmetry_mode == 'GEOM':
            z_sym = (self.l_plies * self.t_thickness) / 2.0
        else:
            z_sym = 0.0
        angle_groups = {}
        for obj in grid.values():
            rp = obj.representative_point().coords[0]
            for i, angle in enumerate(obj.angle):
                for seq in range(self.sequences):
                    seq_offset = len(self.t_angles) * self.t_thickness
                    z_coord = (
                        z_sym + ((2 * i + 1) / 2.0) * self.t_thickness + seq * seq_offset)
                    # print i, seq, angle, [rp[0], rp[1], z_coord]
                    if angle not in angle_groups:
                        angle_groups[angle] = [[rp[0], rp[1], z_coord], ]
                    else:
                        angle_groups[angle].append(
                            [rp[0], rp[1], z_coord])
        if self.symmetry_mode == 'GEOM':
            for key, coord_list in angle_groups.iteritems():
                angle_groups[key].extend(
                    [(x, y, 2.0 * z_sym - z) for (x, y, z) in coord_list])
        return angle_groups

    def assign_properties_3d(self, part, part_grid):
        '''
        Method to assign material properties to the partitioned regions of a 3D
        solid part.
        '''
        cells = part.cells
        angle_groups = self.group_by_angle_3d(part_grid)
        cell_region = regionToolset.Region(cells=cells[:])
        part.SectionAssignment(
            region=cell_region, sectionName='Tape-Elastic',
            offsetType=MIDDLE_SURFACE, thicknessAssignment=FROM_SECTION,
            offset=0.0, offsetField='')

        # assign orientations
        for angle in angle_groups.keys():
            i = 0
            for i, rep_point in enumerate(angle_groups[angle]):
                if i == 0:
                    selected_cells = cells.findAt((rep_point, ))
                else:
                    selected_cells += cells.findAt((rep_point, ))
            region = part.Set(
                cells=selected_cells, name='Cells-{}'.format(int(angle)))
            part.MaterialOrientation(
                region=region, orientationType=SYSTEM, localCsys=None,
                stackDirection=STACK_3, additionalRotationField='',
                angle=angle, additionalRotationType=ROTATION_ANGLE,
                axis=AXIS_3, fieldName='')

    def group_by_layup_2d(self, grid):
        '''
        Method which iterates over polygons to identify the coordinates of 
        areas with same stacking sequence.
        '''
        layup_groups = {}
        for obj in grid.values():
            layup = tuple(obj.angle)
            rp = obj.representative_point().coords[0]
            if layup not in layup_groups:
                layup_groups[layup] = [[rp[0], rp[1], 0.0], ]
            else:
                layup_groups[layup].append([rp[0], rp[1], 0.0])
        return layup_groups

    def assign_properties_2d(self, part, part_grid):
        '''
        Method to assign material properties to a partitioned shell part.
        '''
        faces = part.faces
        layup_groups = self.group_by_layup_2d(part_grid)

        for i, (layup, coord_list) in enumerate(layup_groups.iteritems()):
            for j, rp in enumerate(coord_list):
                if j == 0:
                    selected_faces = faces.findAt((rp, ))
                else:
                    selected_faces += faces.findAt((rp, ))
            region = regionToolset.Region(faces=selected_faces)
            composite_layup = part.CompositeLayup(
                name='CompositeLayup-{}'.format(i), description='',
                elementType=SHELL, offsetType=MIDDLE_SURFACE,
                symmetric=True, thicknessAssignment=FROM_SECTION)
            composite_layup.Section(
                preIntegrate=OFF, integrationRule=SIMPSON, useDensity=OFF,
                thicknessType=UNIFORM, poissonDefinition=DEFAULT,
                temperature=GRADIENT)
            composite_layup.ReferenceOrientation(
                orientationType=GLOBAL, localCsys=None, axis=AXIS_3,
                fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0)
            for sequence in range(self.sequences):
                for n, angle in enumerate(layup):
                    name_of_ply = 'Ply-{}-{}'.format(n, sequence)
                    composite_layup.CompositePly(
                        material='Tape-Elastic', numIntPoints=1, axis=AXIS_3,
                        thicknessType=SPECIFY_THICKNESS, angle=0.0,
                        thickness=self.t_thickness, orientationValue=angle,
                        orientationType=SPECIFY_ORIENT, suppressed=False,
                        additionalRotationField='', plyName=name_of_ply,
                        additionalRotationType=ROTATION_NONE, region=region)

    def mesh_part_3d(self, part):
        '''
        Method to mesh a 3D part using reduced order C3D8R elements with
        enhanced hourglass control and distortion control.
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

    def mesh_part_2d(self, part):
        '''
        Method used to mesh a 2D part using shell elements (with enhanced
        hourglass control).
        '''
        face_region = regionToolset.Region(faces=part.faces[:])
        elemType1 = mesh.ElemType(
            elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF,
            hourglassControl=ENHANCED)
        elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
        part.setElementType(
            regions=face_region, elemTypes=(elemType1, elemType2))
        part.seedPart(size=self.shell_mesh)
        part.generateMesh()


if __name__ == '__main__':
    # Material parameters
    tape_angles = (0, 45, 90, -45)   # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.213  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.104  # ratio of undulation amplitude to length
    number_of_plies = 8  # symmetric 4 ply laminate

    # RVE dimensions
    x_min = 0.0
    x_max = 75.0
    y_min = 0.0
    y_max = 50.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # # Test 3D model creation
    # mdl = AP_PLY_Elastic('TestModel')
    # mdl.symmetry_mode = 'GEOM'
    # mdl.solid_mesh = 1.5  # set mesh size for testing purposes
    # mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
    #                             cured_ply_thickness, undulation_ratio,
    #                             number_of_plies, x_shift=0.0, y_shift=0.0)
    # mdl.create_materials()
    # part = mdl.create_part_3d(specimen_size)
    # grid, lines = mdl.create_grid_elastic(specimen_size)
    # mdl.partition_part(part, lines, cells=True)
    # mdl.partition_through_thickness(part)
    # mdl.place_tows_elastic(grid, specimen_size)
    # mdl.assign_properties_3d(part, grid)
    # mdl.mesh_part_3d(part)

    # Test 2D model creation
    mdl = AP_PLY_Elastic('TestModel')
    mdl.symmetry_mode = 'BC'
    mdl.shell_mesh = 1.5  # set mesh size for testing purposes
    mdl.solid_mesh = 1.5
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies, x_shift=0.0, y_shift=0.0)
    mdl.create_materials()
    part = mdl.create_part_2d(specimen_size)
    grid, lines = mdl.create_grid_elastic(specimen_size)
    mdl.partition_part(part, lines, cells=False)
    mdl.place_tows_elastic(grid, specimen_size)
    mdl.assign_properties_2d(part, grid)
    mdl.mesh_part_2d(part)
