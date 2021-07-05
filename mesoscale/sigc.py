''' 
Interlaced Laminate Analysis: Shapely Interlaced Geometry Creation

This script is used to define the geometry of interlaced laminates.
Geometries are created using the Python Shapely library. The module is
imported by the tape_placement module to create geometries in Abaqus.

(c) Rutger Kok, 23/11/2020
'''
from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from shapely.ops import unary_union, polygonize
from math import cos, radians, tan
from collections import OrderedDict
import pickle
import os
from itertools import groupby
from random import uniform

# initialize additional polygon attributes
setattr(Polygon, 'angle', None)  # create attribute in Polygon class
setattr(Polygon, 'object_type', None)
setattr(Polygon, 'und_pairs', None)


class Interlaced():
    def __init__(self, tape_angles, tape_widths, undulation_width=1.0,
                 specimen=None):
        '''
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
        x, y = zip(*self.specimen.exterior.coords)
        x_max = max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)
        # to ensure the objects on the edges of the specimen are properly
        # determined the specimen size must be increased initially, the excess
        # is later trimmed using the trim_to_specimen function.
        self.specimen_buffered = self.specimen.buffer(2.5, join_style=2)

        # create polygon grids
        # name grid so that it can be saved to reduce future processing time
        grid_name = '{}_{}_{}_{}'.format(
            '-'.join(map(str, self.t_angles)),
            self.t_widths[0],
            self.u_width,
            '-'.join(map(str, (x_min, y_min, x_max, y_max))))
        # define grid storage directory
        root = 'C:\\GitHub\\interlaced_model_creation\\mesoscale\\grids'
        self.grid_path = '{}\\{}.p'.format(root, grid_name)
        # if grid file exists load the file
        if os.path.isfile(self.grid_path):
            print 'Loading grid from file...'
            self._grids = pickle.load(open(self.grid_path, "rb"))
        else:  # if no grid file exists then call the gri creation function
            print 'Creating grid...'
            self._grids = self.create_grids()

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
        x_max = max(x)  # determine max x-value for offset
        y_max = max(y)

        # define mirror point (used to mirror the Polygon boundaries)
        mirror_point = Point([(0.0, 0.0), (0.0, 0.0)])

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
        polygon_dict = {idx: polygon for (idx, polygon) in enumerate(trimmed)}
        # remove sliver polygons
        min_area = 1.0e-4
        clean_dict = {
            idx: polygon for (idx, polygon) in polygon_dict.iteritems()
            if polygon.area > min_area}
        # save the layer grids dictionary for faster run-times in the future
        pickle.dump(clean_dict, open(self.grid_path, "wb"))
        return clean_dict

    def make_pass(self, pass_angle=0, pass_coords=None):
        '''
        Function used to place a single tape/tow of an interlaced laminate.
        In other words, this assigns object types etc. to the Polygon objects
        in self._grids.

        Args:
            pass_angle (float): Angle of the tape/tow to be placed
            pass_coords (tuple): Coordinates of the tape/tow to be placed.
                NOTE: a tape/tow should be specified using its coordinates when
                horizontal and then rotated. Do not rotate outside of the
                Interlaced class

        Returns:
            None
        '''
        # set tape dimensions or use default
        if pass_coords is None:
            w = self.t_widths[0]
            pass_coords = ([(-300.0, (w / 2) - self.u_width),
                            (300.0, (w / 2) - self.u_width),
                            (300.0, (-w / 2) + self.u_width),
                            (-300.0, (-w / 2) + self.u_width)])

        # define tape boundaries
        tol = 1e-6  # tolerance for buffering
        pass_bounds = Polygon(pass_coords)
        pass_rotated_bounds = rotate(
            pass_bounds, pass_angle).buffer(tol, join_style=2)
        # identify which Polygons from the create_grids function are within
        # the boundaries of the rotated tape
        pass_path = [(obj_id, obj) for (obj_id, obj)
                     in self._grids.iteritems()
                     if obj.within(pass_rotated_bounds)]
        if pass_path:
            for obj_id, _ in pass_path:
                self.set_angle(self._grids[obj_id], pass_angle)
                self.set_type(self._grids[obj_id], 'Tape')
            self.create_resin_regions(pass_path, pass_angle)
            return pass_path, pass_angle
        else:
            return [], pass_angle

    def create_resin_regions(self, pass_path, pass_angle):
        '''
        This method creates the "resin regions" which border tape/tow regions.
        For each pass_path, the tape/tow regions in each layer are merged.
        The merged regions are buffered to determine the Polygon objects
        which should be assigned the "Resin" object type.

        Args:
            pass_paths (list): a list of pass_paths (which are lists of
                Polygon objects).

        Returns:
            None
        '''
        tol = 1.0e-6
        # create nested list of objects in pass_path by layer number
        merged_tape = unary_union(list(zip(*pass_path))[1])
        merged_tape_buff = merged_tape.buffer(
            self.u_width + tol, join_style=2)
        search_area = merged_tape_buff.difference(merged_tape)
        resin_regions = [(obj_id, obj) for (obj_id, obj)
                         in self._grids.iteritems()
                         if obj.within(search_area)]
        for (obj_id, obj) in resin_regions:
            self.set_angle(self._grids[obj_id], pass_angle)
            self.set_type(self._grids[obj_id], 'Resin')

    def create_undulations(self, pass_paths, pass_angles):
        '''
        This method deterines the regions in a laminate where tows are
        undulating from one layer to the next.

        All objects in a pass_path are sorted into lists corresponding to
        their layer number. Starting from the bottom (layer 1) we buffer the
        objects in the next layer by twice the undulation width, then check
        which objects in the current layer are within that buffered area.
        These regions are identified as Undulation regions.

        Args:
            pass_paths (list): a list of pass_paths (which are lists of
                Polygon objects).
            pass_angles (list): A list of the angles for the respective
                pass_paths.

        Returns:
            new_paths (list): same as the pass_paths argument but with
                updated object types.
        '''
        tol = 1.0e-6  # tolerance for buffering
        for i, pass_path in enumerate(pass_paths):
            # create nested list of objects in pass_path by layer number
            obj_by_layer = [
                [(obj_id, obj) for (obj_id, obj) in pass_path
                 if obj.angle.index(pass_angles[i]) == s]
                for s in range(len(self.t_angles))]
            search_area = None
            for n in range(1, len(obj_by_layer)):
                # identify objects in current layer and in the layer above
                top_union = unary_union(
                    [obj for (obj_id, obj) in obj_by_layer[n]])
                bottom_union = unary_union(
                    [obj for (obj_id, obj) in obj_by_layer[n - 1]])
                top_union_buff = top_union.buffer(self.u_width, join_style=2)
                bottom_union_buff = bottom_union.buffer(
                    self.u_width, join_style=2)
                search_area = top_union_buff.intersection(bottom_union_buff)
                search_area = search_area.buffer(tol, join_style=2)
                # check search area for the undulations
                detected_undulations = []
                if search_area:
                    detected_undulations = [
                        (obj_id, obj) for (obj_id, obj) in pass_path
                        if obj.within(search_area)]
                # assign properties to undulations
                for obj_id, obj in detected_undulations:
                    self.set_type(
                        self._grids[obj_id], 'Undulation', index=n)
                    self.set_type(
                        self._grids[obj_id], 'Undulation', index=n - 1)

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
            if insert:
                obj.angle.insert(0, angle)
            else:
                obj.angle.append(angle)

    def set_type(self, obj, typ, index=None):
        '''
        Args:
            obj (Shapely Polygon): a Polygon object to which we assign an angle
            angle (float): the angle to assign to the Polygon
            insert (boolean): True if angle must be prepended rather than
                appended to the list.
        '''
        if obj.object_type is None:
            obj.object_type = [typ, ]
        else:
            if index is None:
                obj.object_type.append(typ)
            else:
                try:
                    obj.object_type[index] = typ
                except IndexError:
                    pass

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
        for (obj_id, obj) in self._grids.iteritems():
            if obj.within(self.specimen.buffer(1.0e-6, join_style=2)):
                trimmed_grids[obj_id] = obj
            else:
                continue
        return trimmed_grids

    def object_plot(self):
        '''
        Function to plot the contours of different Polygon types (Tape, Resin,
        Undulation) for debugging of tape_placement, and sigc scripts.
        A different subplot is used for each layer in the laminate.

        Args:
            grid (dict): nested dictionary containing Shapely Polygons
            angles (list): angles in interlaced laminate
            object_type (string): the type (e.g. Tape) of object to plot.
            specimen (Shapely Polygon): the exterior dimensions of the specimen
                which will be the exterior dimensions of the plot.

        Returns:
            None
        '''
        from matplotlib.patches import Polygon as matplotlib_polygon
        import matplotlib.pyplot as plt
        grid = self.trim_to_specimen()
        f1, axes = plt.subplots(1, 1)
        sorted_objects = sorted(grid.values(), key=lambda x: x.angle[0])
        groups = {}
        for key, group in groupby(sorted_objects, key=lambda x: x.angle[0]):
            if key is None:
                key = ['None', ]
            groups[key] = list(group)
        for angle, objects in groups.iteritems():
            r = uniform(0, 1)
            g = uniform(0, 1)
            b = uniform(0, 1)
            for poly in objects:
                vertices = poly.exterior.coords
                patch = matplotlib_polygon(
                    vertices, closed=True, edgecolor='black',
                    facecolor=(r, g, b), fill=True, alpha=0.75,
                    label=angle)
                axes.add_patch(patch)
        axes.autoscale()
        axes.set_xlim(-50.0, 50.0)
        axes.set_ylim(-75.0, 75.0)
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), prop={'size': 6})
        plt.show(f1)

    def test_plot(self, geoms, colours):
        # imports in function definition to avoid issues when importing the
        # script in Abaqus
        from matplotlib import pyplot
        from matplotlib.patches import Polygon as matplotlib_polygon
        '''
        This function is a debugging tool. The function creates plots of the
        geometries passed to it as arguments, in the colours specified.

        Args:
            geoms (list): list of Shapely polygons to be plotted.
            colours (list): list of colours (ordered) to apply to the
                different geoms.

        Returns:
            None
        '''
        fig = pyplot.figure(1)
        ax = fig.add_subplot(111)
        # plot the geometries
        for i, g in enumerate(geoms):
            if g:
                if (g.geom_type == 'GeometryCollection' or
                    g.geom_type == 'MultiPolygon'):
                        for geom in g:
                            vertices = geom.exterior.coords
                            patch = matplotlib_polygon(
                                vertices, closed=True, facecolor=colours[i],
                                fill=True, alpha=0.5)
                            ax.add_patch(patch)
                else:
                    if g.exterior:
                        vertices = g.exterior.coords
                        patch = matplotlib_polygon(
                            vertices, closed=True, facecolor=colours[i],
                            fill=True, alpha=0.5)
                        ax.add_patch(patch)
        # plot the polygon grid
        for (idx, poly) in self._grids.iteritems():
            vertices = poly.exterior.coords
            patch = matplotlib_polygon(
                vertices, closed=True, edgecolor='black', fill=False)
            ax.add_patch(patch)
        ax.autoscale()
        pyplot.show()


if __name__ == '__main__':

    # uw = 1.0
    # w = 10.0
    # tape_angles = (0, 45, 90, -45)
    # tape_widths = (w, ) * len(tape_angles)
    # # initialize laminate
    # interlacedLaminate = Interlaced(
    #     tape_angles, tape_widths, undulation_width=uw)
    # path1, angle1 = interlacedLaminate.make_pass(pass_angle=0)
    # path2, angle2 = interlacedLaminate.make_pass(pass_angle=45)
    # path3, angle3 = interlacedLaminate.make_pass(pass_angle=90)
    # path4, angle4 = interlacedLaminate.make_pass(pass_angle=-45)
    # # crds1 = ([(-300.0, (w / 2) - uw - w), (300.0, (w / 2) - uw - w),
    # #           (300.0, (-w / 2) + uw - w), (-300.0, (-w / 2) + uw - w)])
    # # path5, angle5 = interlacedLaminate.make_pass(
    # #     pass_angle=0, pass_coords=crds1)
    # # crds2 = ([(-300.0, (w / 2) - uw - 2 * w), (300.0, (w / 2) - uw - 2 * w),
    # #           (300.0, (-w / 2) + uw - 2 * w),
    # #           (-300.0, (-w / 2) + uw - 2 * w)])
    # # path6, angle6 = interlacedLaminate.make_pass(
    # #     pass_angle=0, pass_coords=crds2)
    # paths = [path1, path2, path3, path4]  # , path5, path6]
    # angles = [angle1, angle2, angle3, angle4]  # , angle5, angle6]
    # paths = interlacedLaminate.create_undulations(paths, angles)
    # grid = interlacedLaminate.trim_to_specimen()
    # interlacedLaminate.object_plot()

    # uw = 1.0
    # w = 10.0
    # tape_angles = (0, 60, -60)
    # tape_widths = (w, ) * len(tape_angles)
    # # initialize laminate
    # interlacedLaminate = Interlaced(
    #     tape_angles, tape_widths, undulation_width=uw)
    # path1, angle1 = interlacedLaminate.make_pass(pass_angle=0)
    # path2, angle2 = interlacedLaminate.make_pass(pass_angle=60)
    # path3, angle3 = interlacedLaminate.make_pass(pass_angle=-60)
    # # crds1 = ([(-300.0, (w / 2) - uw - w), (300.0, (w / 2) - uw - w),
    # #           (300.0, (-w / 2) + uw - w), (-300.0, (-w / 2) + uw - w)])
    # path4, angle4 = interlacedLaminate.make_pass(
    #     pass_angle=0, pass_coords=crds1)
    # paths = [path1, path2, path3]#, path4]
    # angles = [angle1, angle2, angle3]#, angle4]
    # paths = interlacedLaminate.create_undulations(paths, angles)
    # interlacedLaminate.object_plot()

    uw = 1.0
    w = 10.0
    tape_angles = (0, 90)
    tape_widths = (w, ) * len(tape_angles)
    specimen = Polygon(
        [(-20.0, -50.0), (20.0, -50.0), (20.0, 50.0), (-20.0, 50.0)])
    # initialize laminate
    interlacedLaminate = Interlaced(
        tape_angles, tape_widths, undulation_width=uw, specimen=specimen)
    path1, angle1 = interlacedLaminate.make_pass(pass_angle=0)
    crds1 = ([(-300.0, (w / 2) - uw - 2 * w), (300.0, (w / 2) - uw - 2 * w),
              (300.0, (-w / 2) + uw - 2 * w),
              (-300.0, (-w / 2) + uw - 2 * w)])
    path2, angle2 = interlacedLaminate.make_pass(
        pass_angle=0, pass_coords=crds1)
    crds2 = ([(-300.0, (w / 2) - uw - w), (300.0, (w / 2) - uw - w),
             (300.0, (-w / 2) + uw - w), (-300.0, (-w / 2) + uw - w)])
    path3, angle3, id3 = interlacedLaminate.make_pass(
        pass_angle=0, pass_coords=crds2)
    path4, angle4 = interlacedLaminate.make_pass(pass_angle=90)
    paths = [path1, path4]
    angles = [angle1, angle4]
    interlacedLaminate.create_undulations(paths, angles)
    interlacedLaminate.object_plot()
