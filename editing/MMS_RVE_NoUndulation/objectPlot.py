from itertools import permutations
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import cascaded_union


# Function to plot object types by layer
def objPlot(grid, angles, typ, sample=None):
    import matplotlib.pyplot as plt
    # define specimen boundaries
    if sample is None:
        sample = Polygon(
                         [(-50.0, -75.1), (50.0, -75.1), (50.0, 75.1),
                          (-50.0, 75.1)])
    f1, axes = plt.subplots(2, 2, sharex='col', sharey='row')
    f1.suptitle('{} per Layer'.format(typ))

    angleSets = []
    for r in range(0, 3):
        temp = permutations(angles, r)
        for angleCombinations in temp:
            angleSets.append(list(angleCombinations))

    g = 1
    for i in range(2):
        for j in range(2):
            if g <= len(grid):
                tapesInLayer = [tape for (m, tape)
                                in grid[g].iteritems()
                                if tape.objectType == typ]
                merged = []
                for angleSet in angleSets:
                    sameangle = [tape for tape
                                 in tapesInLayer
                                 if tape.angle == angleSet]
                    merged.append((cascaded_union(sameangle), angleSet))
                for angleSet2, angle in merged:
                    if angleSet2.geom_type == 'Polygon':
                        xi, yi = angleSet2.exterior.xy
                        axes[i, j].plot(xi, yi)
                        axes[i, j].text(xi[0], yi[0], angle)
                        axes[i, j].set_title('Layer: {}'.format(g))
                    elif angleSet2.geom_type == 'MultiPolygon':
                        for t in angleSet2:
                            xi, yi = t.exterior.xy
                            colour = np.random.rand(3,)
                            axes[i, j].text(xi[0], yi[0], angle, color=colour)
                            axes[i, j].plot(xi, yi, color=colour)
                axes[i, j].set_title('Layer: {}'.format(g))
                xb, yb = sample.exterior.xy
                axes[i, j].plot(xb, yb)
                g += 1
            else:
                break

    plt.rcParams['axes.grid'] = True
    plt.show(f1)
