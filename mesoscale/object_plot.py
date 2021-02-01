''' Interlaced Laminate Analysis: Object Plot

This module plots Polygon contours using matplotlib. It is used to debug the
sigc and tape_placement modules.

(c) Rutger Kok, 23/11/2020
'''

from itertools import permutations
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import cascaded_union


def object_plot(grid, angles, object_type, specimen=None):
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
    import matplotlib.pyplot as plt
    # define specimen boundaries
    if specimen is None:
        specimen = Polygon(
                         [(-50.0, -75.1), (50.0, -75.1), (50.0, 75.1),
                          (-50.0, 75.1)])

    f1, axes = plt.subplots(2, 2, sharex='col', sharey='row')
    f1.suptitle('{} per Layer'.format(object_type))
    angle_sets = []
    for r in range(0, 3):
        temp = permutations(angles, r)
        for angle_combinations in temp:
            angle_sets.append(list(angle_combinations))

    g = 1
    for i in range(2):
        for j in range(2):
            if g <= len(grid):
                tapes_in_layer = [tape for (m, tape)
                                  in grid[g].iteritems()
                                  if tape.object_type == object_type]
                merged = []
                for angle_set in angle_sets:
                    same_angle = [tape for tape
                                  in tapes_in_layer
                                  if list(tape.angle) == angle_set]
                    merged.append((cascaded_union(same_angle), angle_set))
                for angle_set_2, angle in merged:
                    if angle_set_2.geom_type == 'Polygon':
                        xi, yi = angle_set_2.exterior.xy
                        axes[i, j].plot(xi, yi)
                        axes[i, j].text(xi[0], yi[0], angle)
                        axes[i, j].set_title('Layer: {}'.format(g))
                    elif angle_set_2.geom_type == 'MultiPolygon':
                        for t in angle_set_2:
                            xi, yi = t.exterior.xy
                            colour = np.random.rand(3,)
                            axes[i, j].text(xi[0], yi[0], angle, color=colour)
                            axes[i, j].plot(xi, yi, color=colour)
                axes[i, j].set_title('Layer: {}'.format(g))
                xb, yb = specimen.exterior.xy
                axes[i, j].plot(xb, yb)
                g += 1
            else:
                break

    plt.rcParams['axes.grid'] = True
    plt.show(f1)
