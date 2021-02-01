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
        paths (list): List of tape/tow connectivity for a tape.
        trimmed_grid (dict): Nested dictionary defining the grid of
            tape/undulation objects.
    """

    # create interlaced laminate instance
    laminate = sigc.Interlaced(tape_angles, tape_widths, undulation_width,
                               specimen=specimen_size)

    x, y = zip(*specimen_size.exterior.coords)
    x_max = max(x)  # determine max x-value for offsets
    y_max = max(y)

    paths = []
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
                    tape_coords = ([(-500.0 + ofs + shift, (w / 2.0) - uw),
                                   (500.0 + ofs + shift, (w / 2.0) - uw),
                                   (500.0 + ofs + shift, (-w / 2.0) + uw),
                                   (-500.0 + ofs + shift, (-w / 2.0) + uw)])
                    created_tape = laminate.make_pass(pass_coords=tape_coords,
                                                      pass_angle=90)
                    paths.append(created_tape)
            else:  # for tape/tow placement angles other than 90
                offset = w / cos(radians(a))
                spaced_offset = ((1 + s) * w) / cos(radians(a))
                max_offset = ((y_max + w / 2.0) * cos(radians(a))
                              +  x_max * tan(radians(a)))
                number_offsets = int(max_offset / spaced_offset)
                offset_list = [spaced_offset * k for k
                               in range(-number_offsets - shift_number - 10,
                                        number_offsets - shift_number + 10)]
                shift = shift_number * offset
                for ofs in offset_list:
                    tape_coords = ([(-500.0, ((w / 2.0) - uw) + ofs + shift),
                                   (500.0, ((w / 2.0) - uw) + ofs + shift),
                                   (500.0, ((-w / 2.0) + uw) + ofs + shift),
                                   (-500.0, ((-w / 2.0) + uw) + ofs + shift)])
                    created_tape = laminate.make_pass(pass_coords=tape_coords,
                                                      pass_angle=a)
                    paths.append(created_tape)
    trimmed_grid = laminate.trim_to_specimen()
    return paths, trimmed_grid


if __name__ == '__main__':
    # This section of the script only runs if this script is run directly
    # (i.e. as long as this script is not imported). It plots the geometries
    # for testing purposes.

    from object_plot import object_plot
    # import matplotlib.pyplot as plt

    # define tape width and thickness
    tape_angles = (0, 90)
    tape_widths = (20,) * len(tape_angles)
    tape_spacing = 1
    cpt = 0.18
    undulation_ratio = 0.18
    rw = cpt / undulation_ratio

    x_min = y_min = -(tape_widths[0] / 2.0)
    x_max = y_max = x_min + (tape_spacing + 1) * (tape_widths[0])
    number_of_layers = len(tape_angles) * 2.0  # symmetric
    laminate_thickness = number_of_layers * cpt
    z_max = laminate_thickness / 2.0
    z_min = -z_max
    rve_polygon = Polygon([(x_min, y_min), (x_max, y_min),
                           (x_max, y_max), (x_min, y_max)])

    tape_paths, grid = laminate_creation(tape_angles, tape_widths,
                                         tape_spacing, rve_polygon,
                                         undulation_width=rw)

    object_plot(grid, tape_angles, 'Tape', sample=rve_polygon)
    object_plot(grid, tape_angles, 'Resin', sample=rve_polygon)
    object_plot(grid, tape_angles, 'Undulation', sample=rve_polygon)
