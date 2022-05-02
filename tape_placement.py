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

from math import cos, tan, radians
import sigc as sigc
from shapely.geometry import Polygon


def laminate_creation(tape_angles, tape_widths, tape_spacing, specimen_size,
                      x_shift=0.0, y_shift=0.0, undulation_width=1.0,
                      debug=False):
    '''
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
            AP-PLY laminate.
        x_shift (float): Shifts AP-PLY architecture in x-direction
        y_shift (float): Shifts AP-PLY architecture in y-direction

    Returns:
        pass_paths (list): List of tape/tow connectivity for a tape.
        trimmed_grid (dict): Nested dictionary defining the grid of
            tape/undulation objects.
    '''

    # create AP-PLY laminate instance
    lm = sigc.AP_PLY(
        tape_angles, tape_widths, undulation_width, x_shift=x_shift,
        y_shift=y_shift, specimen=specimen_size)
    specimen_buffered = specimen_size.buffer(2.5, join_style=2)
    x, y = zip(*specimen_buffered.exterior.coords)
    # NOTE: these max values are different from the self.x_min etc. values
    x_max = max(x)
    x_min = min(x)
    y_max = max(y)
    y_min = min(y)

    x_len = max(abs(lm.x_0 - x_max), abs(lm.x_0 - x_min))
    y_len = max(abs(lm.y_0 - y_max), abs(lm.y_0 - y_min))

    pass_paths = []
    pass_angles = []
    tape_props = zip(tape_angles, tape_widths)
    s = tape_spacing  # reassign tape_spacing for brevity
    uw = undulation_width  # reassign undulation_width for brevity
    for n_shift in range(s + 1):
        for (a, w) in tape_props:
            a_rad = radians(a)
            if a == 90.0 or a == -90.0:
                offset = w
                spaced_offset = (1 + s) * w
                # largest offset within specimen
                max_offset = x_len
            else:
                offset = w / cos(a_rad)
                spaced_offset = ((1 + s) * w) / cos(a_rad)
                max_offset = (y_len + (w / 2.0) + abs(x_len * tan(a_rad)))
            n_offsets = int(max_offset / spaced_offset) + 1
            offset_list = [
                spaced_offset * k for k in range(-n_offsets, n_offsets)]
            shift = n_shift * offset
            for ofs in offset_list:
                if a == 90.0 or a == -90.0:
                    x_1 = lm.x_0 - x_len + ofs + shift
                    x_2 = lm.x_0 + x_len + ofs + shift
                    y_1 = ((w / 2.0) - uw)
                    y_2 = ((-w / 2.0) + uw)
                else:
                    x_1 = lm.x_0 - x_len
                    x_2 = lm.x_0 + x_len
                    y_1 = lm.y_0 + ((w / 2.0) - uw) + ofs + shift
                    y_2 = lm.y_0 + ((-w / 2.0) + uw) + ofs + shift
                tape_coords = ([(x_1, y_1), (x_2, y_1),
                                (x_2, y_2), (x_1, y_2)])
                pass_path, pass_angle = lm.make_pass(
                    pass_coords=tape_coords, pass_angle=a)
                pass_paths.append(pass_path)
                pass_angles.append(pass_angle)
    non_empty_index = [i for i, _ in enumerate(pass_paths) if _ != []]
    pass_paths = [pass_paths[j] for j in non_empty_index]
    pass_angles = [pass_angles[k] for k in non_empty_index]
    lm.create_undulations(pass_paths, pass_angles)
    if debug:
        lm.object_plot()  # only for debugging
    return pass_paths, pass_angles, lm.trim_to_specimen()


if __name__ == '__main__':
    # This section of the script only runs if this script is run directly
    # (i.e. as long as this script is not imported). It plots the geometries
    # for testing purposes.

    # define tape width and thickness
    tape_angles = (0, 45, 90, -45)
    # tape_angles = (0, 90)
    tape_widths = (12.7,) * len(tape_angles)
    tape_spacing = 1
    cpt = 0.213
    undulation_ratio = 0.104
    rw = (cpt / undulation_ratio) / 2.0
    x_ext = 29.0
    y_ext = 29.0
    specimen = Polygon(
        [(0.0, 0.0), (x_ext, 0.0), (x_ext, y_ext), (0.0, y_ext), (0.0, 0.0)])

    shift_x = x_ext / 2.0
    shift_y = y_ext / 2.0
    _, _, grid = laminate_creation(
        tape_angles, tape_widths, tape_spacing, specimen,
        x_shift=-shift_x, y_shift=-shift_y, undulation_width=rw, debug=True)
