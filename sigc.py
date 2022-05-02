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

from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from shapely.ops import unary_union, polygonize
from math import cos, radians, tan, sin
from collections import OrderedDict
from itertools import groupby
from random import uniform


class AP_PLY():
    def __init__(self, tape_angles, tape_widths, undulation_width=1.0,
                 x_shift=0.0, y_shift=0.0, specimen=None):
        '''
        Initializes an instance of an AP-PLY laminate.

        Args:
            tape_angles (tuple): Angles of the tapes/tows in the laminate.
            tape_widths (tuple): Widths of the tapes/tows in the laminate.
            specimen (Shapely Polygon): Polygon defining the dimensions of
                the specimen to be modeled.
            undulation_width (float): Width of "undulation regions" in the
                AP-PLY laminate.
            x_shift (float): Shifts the AP-PLY architecture in the x-direction
            y_shift (float): Shifts the AP-PLY architecture in the y-direction

        Returns:
            None
        '''
        self.t_angles = tape_angles
        self.t_widths = tape_widths
        self.u_width = undulation_width
        self.x_shift = x_shift
        self.y_shift = y_shift
        # define specimen boundaries
        if specimen is None:
            self.specimen = Polygon(
                [(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])
        else:
            self.specimen = specimen
        x, y = zip(*self.specimen.exterior.coords)
        self.x_max = max(x)
        self.x_min = min(x)
        self.y_max = max(y)
        self.y_min = min(y)
        x_mid = (self.x_max + self.x_min) / 2.0
        y_mid = (self.y_max + self.y_min) / 2.0
        # define origin point (where pattern starts)
        self.x_0 = x_mid + self.x_shift
        self.y_0 = y_mid + self.y_shift
        # to ensure the objects on the edges of the specimen are properly
        # determined the specimen size must be increased initially, the excess
        # is later trimmed using the trim_to_specimen function.
        self.specimen_buffered = self.specimen.buffer(
            2.0 * self.u_width, join_style=2)
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
        # NOTE: these max values are different from the self.x_min etc. values
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
                    x_1 = self.x_0 + (w / 2.0) - self.u_width + ofs
                    x_2 = self.x_0 + (w / 2.0) + ofs
                    x_3 = self.x_0 + (w / 2.0) + ofs + self.u_width
                    for x_offset in [x_1, x_2, x_3]:
                        coords = [(x_offset, y_max), (x_offset, y_min)]
                        line = LineString(coords)
                        try:
                            line_trim = line.intersection(
                                self.specimen_buffered)
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
                    y_1 = self.y_0 + (w / 2.0) - self.u_width + ofs
                    y_2 = self.y_0 + (w / 2.0) + ofs
                    y_3 = self.y_0 + (w / 2.0) + ofs + self.u_width
                    for y_offset in [y_1, y_2, y_3]:
                        line = self.create_line(
                            y_offset, a_rad, x_min, x_max, ofs)
                        if line:
                            partition_lines.append(line)

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
        return clean_dict

    def create_line(self, y_offset, a_rad, x_min, x_max, offset):
        y_delta = sin(a_rad) * (y_offset - (self.y_0 + offset))
        x_neg = ((((x_min - self.x_0) + y_delta) / cos(a_rad)) + self.x_0)
        x_pos = ((((x_max - self.x_0) + y_delta) / cos(a_rad)) + self.x_0)
        coords = [(x_neg, y_offset), (x_pos, y_offset)]
        line = rotate(
            LineString(coords), a_rad, Point([self.x_0, self.y_0 + offset]),
            use_radians=True)
        try:
            line_trim = line.intersection(self.specimen_buffered)
            return line_trim
        except NotImplementedError:
            return None

    def make_pass(self, pass_angle=0, pass_coords=None):
        '''
        Function used to place a single tape/tow of an AP-PLY laminate.
        In other words, this assigns object types etc. to the Polygon objects
        in self._grids.

        Args:
            pass_angle (float): Angle of the tape/tow to be placed
            pass_coords (tuple): Coordinates of the tape/tow to be placed.
                NOTE: a tape/tow should be specified using its coordinates when
                horizontal and then rotated. Do not rotate outside of the
                AP-PLY class

        Returns:
            None
        '''
        # set tape dimensions or use default
        if pass_coords is None:
            w = self.t_widths[0]
            pass_coords = ([
                (self.x_min, self.y_0 + (w / 2) - self.u_width),
                (self.x_max, self.y_0 + (w / 2) - self.u_width),
                (self.x_max, self.y_0 + (-w / 2) + self.u_width),
                (self.x_min, self.y_0 + (-w / 2) + self.u_width)])

        # define tape boundaries
        tol = 1e-6  # tolerance for buffering
        pass_bounds = scale(Polygon(pass_coords), xfact=5.0)
        pass_rotated_bounds = rotate(
            pass_bounds, pass_angle).buffer(tol, join_style=2)
        # identify which Polygons from the create_grids function are within
        # the boundaries of the rotated tape
        if pass_rotated_bounds.buffer(self.u_width).intersects(self.specimen):
            # self.test_plot([pass_rotated_bounds], ['#2169CF'])
            pass_path = [(obj_id, obj) for (obj_id, obj)
                         in self._grids.iteritems()
                         if obj.within(pass_rotated_bounds)]
        else:
            return [], pass_angle

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
        tol = 1.0e-2
        # create nested list of objects in pass_path by layer number
        merged_tape = unary_union(list(zip(*pass_path))[1])
        merged_tape_buff = merged_tape.buffer(
            self.u_width + tol, join_style=2)
        search_area = merged_tape_buff.difference(
            merged_tape.buffer(-tol, join_style=2))
        # self.test_plot([search_area, ], ['red', ])
        resin_regions = [(obj_id, obj) for (obj_id, obj)
                         in self._grids.iteritems()
                         if obj.within(search_area)]
        for (obj_id, _) in resin_regions:
            self.set_angle(self._grids[obj_id], pass_angle)
            self.set_type(self._grids[obj_id], 'Resin')
            # pass_path.append((obj_id, self._grids[obj_id]))

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
                # self.test_plot([top_union, bottom_union], ['blue', 'green'])
                # check search area for the undulations
                detected_undulations = []
                if search_area:
                    for obj_id, obj in pass_path:
                        if (obj.within(search_area)
                            and (
                                obj.within(top_union)
                                or obj.within(bottom_union))):
                                    detected_undulations.append((obj_id, obj))
                # assign properties to undulations
                for obj_id, obj in detected_undulations:
                    self.set_type(
                        self._grids[obj_id], 'Undulation', index=n)
                    self.set_type(
                        self._grids[obj_id], 'Undulation', index=n - 1)
                    if obj.within(top_union):
                        self.set_und_pairs(
                            self._grids[obj_id], n, obj.angle[n])
                        self.set_und_pairs(
                            self._grids[obj_id], n - 1,
                            (self._grids[obj_id].angle[n - 1], pass_angles[i]))
                    elif obj.within(bottom_union):
                        try:
                            self.set_und_pairs(
                                self._grids[obj_id], n,
                                (obj.angle[n], pass_angles[i]))
                            self.set_und_pairs(
                                self._grids[obj_id], n - 1, obj.angle[n - 1])
                        except:
                            if obj.within(self.specimen):
                                print obj.angle
                            else:
                                continue

    def set_und_pairs(self, obj, index, angles):
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
        if not hasattr(obj, 'und_pairs'):
            setattr(obj, 'und_pairs', [None] * len(self.t_angles))
        obj.und_pairs[index] = angles

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
        if hasattr(obj, 'angle'):
            if insert:
                obj.angle.insert(0, angle)
            else:
                obj.angle.append(angle)
        else:
            setattr(obj, 'angle', [angle, ])

    def set_type(self, obj, typ, index=None):
        '''
        Args:
            obj (Shapely Polygon): a Polygon object to which we assign an angle
            angle (float): the angle to assign to the Polygon
            insert (boolean): True if angle must be prepended rather than
                appended to the list.
        '''
        if hasattr(obj, 'object_type'):
            if index is None:
                obj.object_type.append(typ)
            else:
                try:
                    obj.object_type[index] = typ
                except IndexError:
                    pass
        else:
            setattr(obj, 'object_type', [typ, ])

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

    def object_plot(self, plot_var='angle'):
        '''
        Function to plot the contours of different Polygon types for debugging
        of tape_placement, and sigc scripts.
        A different subplot is used for each layer in the laminate.

        Args:
            plot_var (string): specifies the variable used to sort the
                plotted polygons.

        Returns:
            None
        '''
        from matplotlib.patches import Polygon as matplotlib_polygon
        import matplotlib.pyplot as plt

        grid = self.trim_to_specimen()
        # define variable to plot
        if plot_var == 'angle':
            objects = [obj for obj in grid.values()
                       if hasattr(obj, 'angle')]
            sorted_objects = sorted(objects, key=lambda x: x.angle)
            grouped_objects = groupby(
                sorted_objects, key=lambda x: x.angle)
        elif plot_var == 'und_pairs':
            objects = [obj for obj in grid.values()
                       if hasattr(obj, 'und_pairs')]
            sorted_objects = sorted(objects, key=lambda x: x.und_pairs)
            grouped_objects = groupby(
                sorted_objects, key=lambda x: x.und_pairs)
        elif plot_var == 'object_type':
            objects = [obj for obj in grid.values()
                       if hasattr(obj, 'object_type')]
            sorted_objects = sorted(objects, key=lambda x: x.object_type)
            grouped_objects = groupby(
                sorted_objects, key=lambda x: x.object_type)
        else:
            print 'Not a valid variable'

        f1, axes = plt.subplots(1, 1)
        groups = {}
        for key, group in grouped_objects:
            if key is None:
                key = ['None', ]
            groups[tuple(key)] = list(group)
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
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), prop={'size': 6})
        plt.show(f1)

    def test_plot(self, geoms, colours):
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
        # imports in function definition to avoid issues when importing the
        # script in Abaqus
        import matplotlib.pyplot as plt
        from descartes import PolygonPatch

        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        # plot the geometries
        for i, g in enumerate(geoms):
            if g:
                if (g.geom_type == 'GeometryCollection' or
                    g.geom_type == 'MultiPolygon'):
                        for geom in g:
                            patch = PolygonPatch(
                                geom.buffer(0), alpha=0.5, ec='none',
                                fc=colours[i])
                            ax.add_patch(patch)
                else:
                    if g.exterior:
                        patch = PolygonPatch(
                            g.buffer(0), alpha=0.5, ec='none', fc=colours[i])
                        ax.add_patch(patch)
        # plot the polygon grid
        for (idx, poly) in self.trim_to_specimen().iteritems():
            patch = PolygonPatch(poly.buffer(0), ec='black', fc='none')
            ax.add_patch(patch)
        ax.autoscale()
        ax.set_aspect('equal')
        plt.savefig('test.svg', format='svg')
        plt.show()


if __name__ == '__main__':
    # Purely to test the sigc script
    uw = 1.0
    w = 10.0
    tape_angles = (0, 45, 90, -45)
    tape_widths = (w, ) * len(tape_angles)
    x_ext = 29.0
    y_ext = 29.0
    specimen = Polygon(
        [(0.0, 0.0), (x_ext, 0.0), (x_ext, y_ext), (0.0, y_ext), (0.0, 0.0)])
    # initialize laminate
    ap_ply_laminate = AP_PLY(
        tape_angles, tape_widths, undulation_width=uw, x_shift=0.0,
        y_shift=0.0, specimen=specimen)
    path1, angle1 = ap_ply_laminate.make_pass(pass_angle=45)
    path2, angle2 = ap_ply_laminate.make_pass(pass_angle=-45)
    paths = [path1, path2]
    angles = [angle1, angle2]
    ap_ply_laminate.create_undulations(paths, angles)
    ap_ply_laminate.object_plot()
