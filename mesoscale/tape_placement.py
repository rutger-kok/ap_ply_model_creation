# -*- coding: utf-8 -*-
'''Interlaced Laminate Analysis: Tape Placement

This module defines the geometry of an idealized interlaced laminate.
It is used by interlaced_model which imports this module from within Abaqus.

This module imports the sigc (shapely interlced geometry creation) module to
create a grid of polygons with unassigned object types, then assigns the
object types according to the interlaced lainate parameters passed as
arguments to the laminate_creation function.

(c) Rutger Kok, 23/11/2020
'''
from math import cos, tan, radians
import sigc as sigc
from shapely.geometry import Polygon


def laminate_creation(tape_angles, tape_widths, tape_spacing, specimen_size,
                      undulation_width=1.0):
    """laminate_creation function

    Creates a laminate instance and assigns object types to the laminate based
    on the interlacing parameters.

    Args:
        tape_angles (tuple): Angles of the tapes/tows in the laminate.
        tape_widths (tuple): Widths of the tapes/tows in the laminate.
        tape_spacing (int): Number of tape/tow width gaps between tows placed
            in a single pass.
        specimen_size (Shapely Polygon): Polygon defining the dimensions of
            the specimen to be modeled.
        undulation_width (float): Width of "undulation regions" in the
            interlaced laminate.

    Returns:
        pass_paths (list): List of tape/tow connectivity for a tape.
        trimmed_grid (dict): Nested dictionary defining the grid of
            tape/undulation objects.
    """

    # create interlaced laminate instance
    laminate = sigc.Interlaced(tape_angles, tape_widths, undulation_width,
                               specimen=specimen_size)
    specimen_buffered = specimen_size.buffer(2.5, join_style=2)
    x, y = zip(*specimen_buffered.exterior.coords)
    x_max = max(x) + 2.0 * max(tape_widths)  # determine max x-value for offset
    y_max = max(y) + 2.0 * max(tape_widths)
    d_max = (x_max**2.0 + y_max**2.0)**0.5

    pass_paths = []
    pass_angles = []
    tape_props = zip(tape_angles, tape_widths)
    s = tape_spacing  # reassign tape_spacing for brevity
    uw = undulation_width  # reassign undulation_width for brevity
    for shift_number in range(s + 1):
        for (a, w) in tape_props:
            if a == 90:  # if placement angle is 90 degrees
                offset = w
                # offset equal to width of tape/tow multipled by spacing
                spaced_offset = (1 + s) * offset
                max_offset = x_max + w / 2.0  # largest offset within specimen
                number_offsets = int(max_offset / spaced_offset)
                offset_list = [spaced_offset * k for k
                               in range(-number_offsets - shift_number - 10,
                                        number_offsets - shift_number + 10)]
                shift = shift_number * offset  # offset from starting point
                for ofs in offset_list:
                    tape_coords = ([(-d_max + ofs + shift, (w / 2.0) - uw),
                                   (d_max + ofs + shift, (w / 2.0) - uw),
                                   (d_max + ofs + shift, (-w / 2.0) + uw),
                                   (-d_max + ofs + shift, (-w / 2.0) + uw)])
                    pass_path, pass_angle = laminate.make_pass(
                        pass_coords=tape_coords, pass_angle=a)
                    pass_paths.append(pass_path)
                    pass_angles.append(pass_angle)
            else:  # for tape/tow placement angles other than 90
                offset = w / cos(radians(a))
                spaced_offset = ((1 + s) * w) / cos(radians(a))
                max_offset = ((y_max + w / 2.0) * cos(radians(a))
                              + x_max * tan(radians(a)))
                number_offsets = int(max_offset / spaced_offset)
                offset_list = [spaced_offset * k for k
                               in range(-number_offsets - shift_number - 10,
                                        number_offsets - shift_number + 10)]
                shift = shift_number * offset
                for ofs in offset_list:
                    tape_coords = ([(-d_max, ((w / 2.0) - uw) + ofs + shift),
                                   (d_max, ((w / 2.0) - uw) + ofs + shift),
                                   (d_max, ((-w / 2.0) + uw) + ofs + shift),
                                   (-d_max, ((-w / 2.0) + uw) + ofs + shift)])
                    pass_path, pass_angle = laminate.make_pass(
                        pass_coords=tape_coords, pass_angle=a)
                    pass_paths.append(pass_path)
                    pass_angles.append(pass_angle)
    non_empty_index = [i for i, _ in enumerate(pass_paths) if _ != []]
    pass_paths = [pass_paths[j] for j in non_empty_index]
    pass_angles = [pass_angles[k] for k in non_empty_index]
    laminate.create_undulations(pass_paths, pass_angles)
    # laminate.object_plot()
    return pass_paths, pass_angles, laminate.trim_to_specimen()


if __name__ == '__main__':
    # This section of the script only runs if this script is run directly
    # (i.e. as long as this script is not imported). It plots the geometries
    # for testing purposes.

    # define tape width and thickness
    tape_angles = (0, 45, 90, -45)
    tape_widths = (10.0,) * len(tape_angles)
    tape_spacing = 3
    cpt = 0.18
    undulation_ratio = 0.09
    rw = (cpt / undulation_ratio) / 2.0
    specimen = Polygon(
        [(-20.0, -50.0), (20.0, -50.0), (20.0, 50.0), (-20.0, 50.0)])
    tape_paths, tape_angles, grid = laminate_creation(
        tape_angles, tape_widths, tape_spacing, specimen,
        undulation_width=rw)

    # # define tape width and thickness
    # tape_angles = (0, 90)
    # tape_widths = (10.0,) * len(tape_angles)
    # tape_spacing = 3
    # cpt = 0.18
    # undulation_ratio = 0.09
    # rw = (cpt / undulation_ratio) / 2.0
    # specimen = Polygon(
    #     [(-20.0, -50.0), (20.0, -50.0), (20.0, 50.0), (-20.0, 50.0)])
    # tape_paths, tape_angles, grid = laminate_creation(
    #     tape_angles, tape_widths, tape_spacing, specimen, undulation_width=rw)
