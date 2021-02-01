''' Interlaced Laminate Analysis: Shapely Interlaced Geometry Creation

This script is used to define the geometry of interlaced laminates.
Geometries are created using the Python Shapely library. The module is
imported by the tapePlacement module to create geometries in Abaqus.

(c) Rutger Kok, 23/11/2020
'''
from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from shapely.ops import cascaded_union, unary_union, polygonize
from copy import deepcopy
from math import cos, radians, tan
from itertools import combinations

# initialize additional polygon attributes
setattr(Polygon, 'object_type', None)
setattr(Polygon, 'angle', None)  # create attribute in Polygon class
setattr(Polygon, 'layer', 1)  # create attribute in Polygon class
setattr(Polygon, 'object_id', None)


class Interlaced():
    def __init__(self, tape_angles, tape_widths, undulation_width=1.0,
                 specimen=None):
        ''' Interlaced.__init__ method

        Initializes an instance of an interlaced laminate.

        Args:
            tape_angles (tuple): Angles of the tapes/tows in the laminate.
            tape_widths (tuple): Widths of the tapes/tows in the laminate.
            specimen (Shapely Polygon): Polygon defining the dimensions of
                the specimen to be modeled.
            undulation_width (float): Width of "undulation regions" in the
                interlaced laminate.

        Returns:
            None
        '''
        self.t_angles = tape_angles
        self.t_widths = tape_widths
        self.u_width = undulation_width

        # define specimen boundaries
        if specimen is None:
            self.specimen = Polygon(
                [(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])
        else:
            self.specimen = specimen
        self.specimen_buffered = self.specimen.buffer(2.5, join_style=2)

        # create Polygon grids
        self._grids = self.create_grids()

    def make_pass(self, pass_angle=0, pass_coords=None):
        '''
        Function used to place a single tape/tow of an interlacd laminate.
        In other words, this assigns object types etc. to the Polygon objects
        in self._grids.

        Args:
            pass_angle (float): Angle of the tape/tow to be placed
            pass_coords (tuple): Coordinates of the tape/tow to be placed.
                NOTE: a tape/tow should be specified using its coordinates when
                horizontal and then rotated. Do not rotate outside of the
                Interlaced class

        Returns:
            connectivity (list): A list containing the layer numbers and object
                IDs of all the objects in the pass_path (used to tie them
                together in Abaqus)
        '''

        # set tape dimensions or use default
        if pass_coords is None:
            w = 30.0
            pass_coords = ([(-300.0, (w / 2) - self.u_width),
                            (300.0, (w / 2) - self.u_width),
                            (300.0, (-w / 2) + self.u_width),
                            (-300.0, (-w / 2) + self.u_width)])

        # define tape boundaries
        pass_bounds = Polygon(pass_coords)
        pass_rotated_bounds = rotate(
            pass_bounds, pass_angle).buffer(1 * 10**-8, join_style=2)
        # identify which Polygons from the create_grids function are within
        # the boundaries of the rotated tape
        pass_path = [(obj_id, obj) for (obj_id, obj)
                     in self._grids[1].iteritems()
                     if obj.within(pass_rotated_bounds)]

        pass_path = self.raise_tapes(pass_path, pass_angle)
        pass_path = self.create_undulations(pass_path, pass_angle)
        self.create_resin_regions(pass_path, pass_angle)
        self.grids = self.trim_to_specimen()
        connectivity = self.define_connectivity(pass_path)
        return connectivity

    def create_grids(self):
        '''
        This method splits the Polygon which defines the dimensions of the
        specimen (self.specimen_buffered) into a number of smaller polygons,
        which are bounded by the dimensions of the tapes/tows to be placed.

        Args:
            None

        Returns:
            layer_grids (dict): a nested dictionary containing Polygon
                objects and their IDs for each layer in the laminate.
        '''

        x, y = zip(*self.specimen_buffered.exterior.coords)
        x_max = max(x)  # determine max x-value for offsets
        y_max = max(y)

        # define mirror point (used to mirror the Polygon boundaries)
        mirror_point = Point([(0.0, 0.0), (0.0, 0.0)])

        number_layers = len(self.t_angles)  # number of layers in the laminate
        tapes = zip(self.t_angles, self.t_widths)
        partition_lines = []
        for (a, w) in tapes:
            if a == 90:
                offset = w
                max_offset = x_max - w / 2.0
                number_offsets = int(max_offset / offset)
                offset_list = [offset * k for k in range(number_offsets + 1)]
                for ofs in offset_list:
                    # tape/tow boundaries
                    coords_1 = [((w / 2.0) - self.u_width + ofs, 300.0),
                                ((w / 2.0) - self.u_width + ofs, -300.0)]
                    # boundaries of undulation/resin regions
                    coords_2 = [((w / 2.0) + ofs, 300.0),
                                ((w / 2.0) + ofs, -300.0)]
                    coords_3 = [((w / 2.0) + ofs + self.u_width, 300.0),
                                ((w / 2.0) + ofs + self.u_width, -300.0)]
                    line_1 = LineString(coords_1)
                    line_2 = LineString(coords_2)
                    line_3 = LineString(coords_3)
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
                    coords_1 = [(-300.0, (w / 2.0) - self.u_width + ofs),
                                (300.0, (w / 2.0) - self.u_width + ofs)]
                    coords_2 = [(-300.0, (w / 2.0) + ofs),
                                (300.0, (w / 2.0) + ofs)]
                    coords_3 = [(-300.0, (w / 2.0) + ofs + self.u_width),
                                (300.0, (w / 2.0) + ofs + self.u_width)]
                    rotation_point = Point([(0.0, ofs), (0.0, ofs)])
                    line_1 = rotate(LineString(coords_1), a, rotation_point)
                    line_2 = rotate(LineString(coords_2), a, rotation_point)
                    line_3 = rotate(LineString(coords_3), a, rotation_point)
                    partition_lines.extend([line_1, line_2, line_3])
                reflected_lines = [rotate(line, 180.0, origin=mirror_point)
                                   for line in partition_lines]
                partition_lines = partition_lines + reflected_lines

        # collection of individual linestrings for splitting in a list and add
        # the polygon lines to it.
        partition_lines.append(self.specimen_buffered.boundary)
        partition_lines.append(self.specimen.boundary)
        border_lines = unary_union(partition_lines)
        decomposition = polygonize(border_lines)
        trimmed = [polygon for polygon in decomposition
                   if polygon.within(self.specimen_buffered)]
        polygon_dictionary = {idx: polygon for (idx, polygon)
                              in enumerate(trimmed)}
        layer_grids = {l: deepcopy(polygon_dictionary)
                       for l in range(1, number_layers + 1)}
        return layer_grids

    def trim_to_specimen(self):
        '''
        This method removes polygons outwith the specimen dimensions
        (which are a result of buffering) from the grids dictionary.

        Args:
            None

        Returns:
            trimmed_grids (dict): a nested dictionary containing Polygon
                objects and their IDs for each layer in the laminate.
        '''
        trimmed_grids = {}
        for layer_number, grid_layer in self._grids.iteritems():
            trimmed_grids[layer_number] = {}
            for (obj_id, obj) in grid_layer.iteritems():
                if obj.within(self.specimen):
                    trimmed_grids[layer_number][obj_id] = obj
                else:
                    continue
        return trimmed_grids

    def create_resin_regions(self, pass_path, pass_angle):
        '''
        This method creates the "resin regions" which border tape/tow regions.
        All tape/tow regions in a layer are merged, then these merged regions
        are buffered to determine the Polygon objects which should be assigned
        the "Resin" object type.

        If the regions have already been identified as "Undulation" regions,
        then the angle of the current pass is added to the undulation object
        angle set.

        Args:
            pass_path (list): a list of the tape/tow objects in the current
                pass.
            pass_angle (float): Angle of the tape/tow in the current pass.

        Returns:
            None
        '''
        # create nested list of objects in pass_path by layer number
        obj_by_layer = [[obj for (obj_id, obj) in pass_path if obj.layer == s]
                        for s in range(1, len(self._grids) + 1)]

        for n in range(1, len(obj_by_layer) + 1):
            # merge all objects in a layer
            merged_obj = cascaded_union(obj_by_layer[n - 1])
            # buffer by width of undulations
            merged_obj_buffer = merged_obj.buffer(self.u_width, join_style=2)
            # create list of Polygons within buffered area
            resin_regions = [(obj_id, obj) for (obj_id, obj)
                             in self._grids[n].iteritems()
                             if obj.within(merged_obj_buffer)]
            # assign object type for Polygons in resin_regions
            for (i, (obj_id, obj)) in enumerate(resin_regions):
                if obj.object_type is None:
                    setattr(self._grids[n][obj_id], 'object_type', 'Resin')
                    self._grids[n][obj_id].layer = n
                    self.set_angle(self._grids[n][obj_id], pass_angle)
                    resin_regions[i] = (obj_id, self._grids[n][obj_id])
                elif obj.object_type == 'Undulation' and pass_angle not in obj.angle:
                    self.set_angle(self._grids[n][obj_id], pass_angle,
                                   insert=True)

    def raise_tapes(self, pass_path, pass_angle):
        '''
        This method iterates over each Polygon identified as being
        within the bounds of the current pass (pass_path). The loop checks
        each Polygon's attribute type. If it is 'Polygon' that means that the
        object has not yet been identified as a Tape or Undulation region.
        In this case, the loop sets the object_type attribute of the Polygon
        to 'Tape'.

        If the object_type of a Polygon is 'Tape' that means a previous pass
        has already assigned the 'Tape' object_type to the Polygon in question.
        In other words, there is already a tape placed in this area in this
        layer. If this is the case, the loop sets the attribute of the same
        Polygon in the next grid layer to 'Tape'. Essentially, if a tape
        crossing is detected then the Tape is placed in the next layer.

        The method then calls the in_between method to raise the regions
        between the raised tape regions.

        Args:
            pass_path (list): a list of the tape/tow objects in the current
                pass.
            pass_angle (float): Angle of the tape/tow in the current pass.

        Returns:
            pass_path (list): same as the pass_path method argument but with
                updated object types.
        '''
        for i in range(len(pass_path)):
            obj_id, obj = pass_path[i]
            if obj.object_type is None:
                # if object_type not yet assigned -> assign Tape
                setattr(obj, 'object_type', 'Tape')
                # set angle attibute of Polygon object
                self.set_angle(obj, pass_angle)
            elif obj.object_type == 'Tape':
                # if object_type has already been assigned -> check layer above
                next_layer = 2
                while next_layer <= len(self._grids):
                    next_layer_grid = self._grids[next_layer]
                    if next_layer_grid[obj_id].object_type is None:
                        setattr(next_layer_grid[obj_id], 'object_type', 'Tape')
                        next_layer_grid[obj_id].layer = next_layer
                        self.set_angle(next_layer_grid[obj_id], pass_angle)
                        pass_path[i] = (obj_id, next_layer_grid[obj_id])
                        break
                    else:
                        next_layer += 1
        # raise regions in-between raised tapes
        pass_path = self.in_between(pass_path, pass_angle, 'Tape')
        return pass_path

    def in_between(self, pass_path, pass_angle, objType):
        '''
        This method is called by the raise_tapes method. All tape/tow regions
        in a layer are merged, and then buffered. The intersection of the
        buffered regions identifies the Polygons between raised Tapes which
        need to be raised to the same level.

        In essence it converts a tape profile from this:
        ____---_---___
        to this:
        ____-------___

        Args:
            pass_path (list): a list of the tape/tow objects in the current
                pass.
            pass_angle (float): Angle of the tape/tow in the current pass.

        Returns:
            pass_path (list): same as the pass_path method argument but with
                updated object types.
        '''
        # create nested list of objects in pass_path by layer number
        obj_by_layer = [[obj for (obj_id, obj) in pass_path
                         if obj.layer == s and obj.object_type == objType]
                        for s in range(1, len(self._grids) + 1)]

        for n in range(1, len(obj_by_layer)):
            # raise regions between objects in the same layer
            merged_obj = cascaded_union(obj_by_layer[n])
            try:
                merged_obj_buffer = [obj.buffer(2 * self.u_width, join_style=2)
                                     for obj in merged_obj.geoms]
                merged_obj_convex_difference = cascaded_union(
                    [obj_1.intersection(obj_2) for obj_1, obj_2
                     in combinations(merged_obj_buffer, 2)])
                for (obj_id, obj) in pass_path:
                    if obj.within(merged_obj_convex_difference):
                        obj = self._grids[n + 1][obj_id]
                        obj.object_type = objType
                        obj.layer = n + 1
                        self.set_angle(obj, pass_angle)
                        # idx = pass_path.index(obj_id)
                        idx = [i for (i, el) in pass_path].index(obj_id)
                        pass_path[idx] = (obj_id, obj)
            except AttributeError:
                pass
        return pass_path

    def create_undulations(self, pass_path, pass_angle):
        '''
        This method deterines the regions in a laminate where tows are
        undulating from one layer to the next.

        All objects are sorted into lists corresponding to their layer number.
        Starting from the bottom (layer 1) we buffer the objects in the next
        layer by twice the undulation width, then check which objects in the
        current layer are within that buffered area. These regions are
        identified as Undulation regions.

        Args:
            pass_path (list): a list of the tape/tow objects in the current
                pass.
            pass_angle (float): Angle of the tape/tow in the current pass.

        Returns:
            pass_path (list): same as the pass_path method argument but with
                updated object types.
        '''
        # create nested list of objects in pass_path by layer number
        obj_by_layer = [[(obj_id, obj) for (obj_id, obj)
                         in pass_path if obj.layer == s]
                        for s in range(1, len(self._grids) + 1)]

        for n in range(1, len(obj_by_layer)):
            # merge objects in layer n and buffer
            top_layer = cascaded_union(
                [obj for (obj_id, obj) in obj_by_layer[n]])
            top_layer_buffer = top_layer.buffer(
                2 * self.u_width + 1.0, join_style=2)
            bottom_layer = obj_by_layer[n - 1]
            detected_undulations = [
                (obj_id, obj) for (obj_id, obj) in bottom_layer
                if obj.within(top_layer_buffer)]
            for (obj_id, obj) in detected_undulations:
                bottom_obj = self._grids[n][obj_id]
                top_obj = self._grids[n + 1][obj_id]
                top_obj.object_type = 'Undulation'
                top_obj.layer = n + 1
                bottom_obj.object_type = 'Undulation'
                if len(bottom_obj.angle) != 2:
                    self.set_angle(bottom_obj, pass_angle)
                self.set_angle(top_obj, pass_angle)
                pass_path.append([obj_id, top_obj])

        return pass_path

    def define_connectivity(self, pass_path):
        '''
        This method is used to define the connectivity of tapes/tows in a
        single pass.

        Args:
            pass_path (list): a list of the tape/tow objects in the current
                pass.

        Returns:
            connectivity (list): A list containing the layer numbers and object
                IDs of all the objects in the pass_path (used to connect tie
                them together in Abaqus)
        '''
        connectivity = []
        for (obj_id, obj) in pass_path:
            if obj_id in self.grids[int(obj.layer)].keys():
                ident = '{}-{}'.format(obj.layer, obj_id)
                setattr(obj, 'object_id', ident)
                connectivity.append(ident)
        return connectivity

    def set_angle(self, obj, angle, insert=False):
        '''
        Function to set the angle attribute of Polygon instances.
        Cannot setattr(Polygon, 'angle, set()) because sets are  mutable.
        This would result in all the Polygons sharing the same angle attribute.
        Instead initialize the angle attribute with None and then use set_angle
        to either create a set with the first angle or add to the set if an
        angle has already been assigned.

        Args:
            obj (Shapely Polygon): a Polygon object to which we assign an angle
            angle (float): the angle to assign to the Polygon
            insert (boolean): True if angle must be prepended rather than
                appended to the list.
        '''
        if obj.angle is None:
            obj.angle = [angle, ]
        else:
            if angle not in obj.angle:
                if insert:
                    obj.angle.insert(0, angle)
                else:
                    obj.angle.append(angle)


# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.


if __name__ == '__main__':
    from object_plot import object_plot
    # import matplotlib.pyplot as plt

    tape_angles = (0, 90)
    tape_widths = (30.0, ) * len(tape_angles)

    # initialize laminate
    interlacedLaminate = Interlaced(tape_angles, tape_widths)

    interlacedLaminate.make_pass(pass_angle=0)
    w = 30.0
    uw = 1.0
    crds = ([(-300.0, (w / 2) - uw - w), (300.0, (w / 2) - uw - w),
             (300.0, (-w / 2) + uw - w), (-300.0, (-w / 2) + uw - w)])
    interlacedLaminate.make_pass(pass_angle=0, pass_coords=crds)
    interlacedLaminate.make_pass(pass_angle=90)
    crds = ([(-300.0, (w / 2) - uw - 2 * w), (300.0, (w / 2) - uw - 2 * w),
             (300.0, (-w / 2) + uw - 2 * w), (-300.0, (-w / 2) + uw - 2 * w)])
    interlacedLaminate.make_pass(pass_angle=0, pass_coords=crds)
    grid = interlacedLaminate.trim_to_specimen()

    # for (idx, pol) in grid[1].iteritems():
    #     x, y = pol.exterior.xy
    #     plt.plot(x, y)
    # plt.show()

    object_plot(grid, tape_angles, 'Tape')
    object_plot(grid, tape_angles, 'Resin')
    object_plot(grid, tape_angles, 'Undulation')
