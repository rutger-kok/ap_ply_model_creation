''' Interlaced Laminate Analysis: Object Plot

This module plots Polygon contours using matplotlib. It is used to debug the
sigc and tape_placement modules.

(c) Rutger Kok, 23/11/2020
'''

from itertools import groupby
from shapely.geometry import Polygon
from matplotlib.patches import Polygon as matplotlibPolygon
import matplotlib.pyplot as plt


def object_plot(grid, specimen=None):
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
    # define specimen boundaries
    if specimen is None:
        specimen = Polygon(
            [(-50.0, -75.1), (50.0, -75.1), (50.0, 75.1), (-50.0, 75.1)])

    f1, axes = plt.subplots(2, 2, sharex='col', sharey='row')
    colors = {'Undulation': 'red', 'Resin': 'blue', 'Tape': 'green'}
    g = 1
    for i in range(2):
        for j in range(2):
            if g <= len(grid):
                sorted_objects = sorted(
                    grid[g].values(), key=lambda x: x.object_type)
                groups = {key: list(v for v in valuesiter)
                          for key, valuesiter
                          in groupby(
                              sorted_objects, key=lambda x: x.object_type)}
                for object_type, objects in groups.iteritems():
                    for poly in objects:
                        vertices = poly.exterior.coords
                        if object_type is None:
                            patch = matplotlibPolygon(
                                vertices, closed=True, edgecolor='black',
                                fill=False)
                        else:
                            patch = matplotlibPolygon(
                                vertices, closed=True, edgecolor='black',
                                facecolor=colors[object_type], fill=True,
                                alpha=0.75)
                        axes[i, j].add_patch(patch)
                axes[i, j].autoscale()
                axes[i, j].set_title('Layer: {}'.format(g))
                g += 1
            else:
                break
    plt.show(f1)