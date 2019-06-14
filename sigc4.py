from shapely.geometry import Polygon, Point, LineString, JOIN_STYLE
from shapely import affinity
from shapely.ops import cascaded_union, split, linemerge, unary_union, polygonize
from copy import deepcopy
import numpy as np
import math

# ----------------------------------------------------------------------------
# Function to plot object types by layer

def objPlot(grid, angles, typ, sample=None):
    import matplotlib.pyplot as plt
    # define specimen boundaries
    if sample == None:
        sample = Polygon(
            [(-50.0, -75.1), (50.0, -75.1), (50.0, 75.1), (-50.0, 75.1)])
    f1, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    f1.suptitle('{} per Layer'.format(typ))

    g = 0
    for i in range(3):
        for j in range(3):
            if g <= len(grid)-1:
                tapesInLayer = [tape for (m,tape) 
                        in grid[g].iteritems()
                        if tape.objectType == typ]
                merged = []
                for ang in angles:
                    sameangle = [tape for tape 
                        in tapesInLayer
                        if tape.angle == ang]
                    merged.append(cascaded_union(sameangle))
                for angleSet in merged:           
                    if angleSet.geom_type == 'Polygon':
                        xi,yi = angleSet.exterior.xy
                        axes[i,j].plot(xi,yi)
                        axes[i,j].set_title('Layer: {}'.format(g))
                    elif angleSet.geom_type == 'MultiPolygon':
                        for t in angleSet:
                            xi,yi = t.exterior.xy
                            axes[i,j].plot(xi,yi)
                axes[i,j].set_title('Layer: {}'.format(g))
                xb, yb = sample.exterior.xy
                axes[i,j].plot(xb,yb)
                g += 1
            else:
                break
    # g = 0
    # for i in range(3):
    #     for j in range(3):
    #         if g <= len(grid)-1:
    #             tapesInLayer = [tape for (m,tape) 
    #                     in grid[g].iteritems()
    #                     if tape.objectType == typ]
    #             # if tapesInLayer.geom_type == 'Polygon':
    #             #     xi,yi = tapesInLayer.exterior.xy
    #             #     axes[i,j].plot(xi,yi)
    #             #     axes[i,j].set_title('Layer: {}'.format(g))
    #             # elif tapesInLayer.geom_type == 'MultiPolygon':
    #             for t in tapesInLayer:
    #                 xi,yi = t.exterior.xy
    #                 axes[i,j].plot(xi,yi)
    #             axes[i,j].set_title('Layer: {}'.format(g))
    #             xb, yb = sample.exterior.xy
    #             axes[i,j].plot(xb,yb)
    #             g += 1
    #         else:
    #             break
    # plt.show(f1)
    plt.show(f1)

# ----------------------------------------------------------------------------

# initialize additional polygon attributes
setattr(Polygon, 'objectType', 'Polygon')
setattr(Polygon, 'angle', set()) # create attribute in Polygon class
setattr(Polygon, 'layer', 1) # create attribute in Polygon class

# Function to partition the specimen into a grid of polygons depending on the 
# tape angles of the laminate
def createGrids(tapeAngles, tapeWidths, r=1.0, sample=None):
    # define specimen boundaries
    if sample == None:
        sample = Polygon(
            [(-50.0, -75.1), (50.0, -75.1), (50.0, 75.1), (-50.0, 75.1)])
    
    tapes = zip(tapeAngles, tapeWidths)
    partitionLines = []
    for tape in tapes:
        a = tape[0]  # angle
        w = tape[1]  # width
        if a == 90:
            offset = (w+2*r)
            maxOffset = 50.0-w/2.0
            numberOffsets = int(maxOffset/offset)
            offsetList = [offset*k for k 
                            in range(-numberOffsets-1,numberOffsets+1)]
            for os in offsetList:
                tapeLineCoords = [((w/2.0)+os, 100.0), ((w/2.0)+os, -100.0)]
                resinLineCoords1 = (
                        [((w/2.0)+os+r, 100.0),
                        ((w/2.0)+os+r, -100.0)])
                resinLineCoords2 = (
                        [((w/2.0)+os+2*r, 100.0),
                        ((w/2.0)+os+2*r, -100.0)])

                tapeLine = LineString(tapeLineCoords)
                resinLine1 = LineString(resinLineCoords1)
                resinLine2 = LineString(resinLineCoords2)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])

        else:
            offset = (w+2*r)/math.cos(math.radians(a))
            maxOffset = (75.0 - (w/2.0)*math.cos(math.radians(a))
                    + 50.0*math.tan(math.radians(a)))
            numberOffsets = int(maxOffset/offset)
            offsetList = [offset*k for k 
                            in range(-numberOffsets-1,numberOffsets+1)]
            for os in offsetList:
                tapeLineCoords = [(-100.0, (w/2.0)+os), (100.0, (w/2.0)+os)]
                resinLineCoords1 = (
                        [(-100.0, (w/2.0)+os+r), (100.0, (w/2.0)+os+r)])
                resinLineCoords2 = (
                        [(-100.0, (w/2.0)+os+2*r), (100.0, (w/2.0)+os+2*r)])
                rotPoint = Point([(0.0, os), (0.0, os)])
                tapeLine = affinity.rotate(
                    LineString(tapeLineCoords), a, rotPoint)
                resinLine1 = affinity.rotate(
                    LineString(resinLineCoords1), a, rotPoint)
                resinLine2 = affinity.rotate(
                    LineString(resinLineCoords2), a, rotPoint)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])

    # collection of individual linestrings for splitting in a list and add 
    # the polygon lines to it.
    partitionLines.append(sample.boundary) 
    merged_lines = linemerge(partitionLines)
    border_lines = unary_union(merged_lines)
    decomposition = polygonize(border_lines)     
    trimmed = [plygn for plygn in decomposition if plygn.within(sample)]
    polyDict = {key: value for (key, value) in enumerate(trimmed)}
    layerGrids = [deepcopy(polyDict) for l in range(len(tapeAngles))]
    return layerGrids

# for poly in grids[0]:
#     x,y = poly.exterior.xy
#     plt.plot(x,y)
#     plt.grid(True)
# plt.show()

# Class defining a single pass of a an ATP machine
# inputs: grid of the specimen, angle of the pass, coords, width of resin
# NOTE: a pass should be specified using its coordinates when horizontal 
# and then rotated. Do not rotate outside of the machinePass class
class machinePass():
    def __init__(self, grids, lNum=1, angle=0, coords=None, resinW=1.0):
        self.layer = lNum
        self.angle = angle
        if coords == None:
            w = 10.0
            self.coords = (
                    [(-100.0, w/2), (100.0, w/2), (100.0, -w/2), (-100.0, -w/2)]
                    )
        else:
            self.coords = coords
        
        bounds = Polygon(self.coords)
        rotatedBounds = affinity.rotate(
                bounds, self.angle).buffer(1*10**-8, join_style=2)

        self.tapePath = [(i,plygn) for (i,plygn) 
                            in grids[self.layer-1].iteritems() 
                            if plygn.within(rotatedBounds)]

        # lift up overlapping regions
        for index in range(len(self.tapePath)):
            i, obj = self.tapePath[index]
            if obj.objectType == 'Polygon':
                obj.objectType = 'Tape'
                obj.angle = self.angle
            elif obj.objectType == 'Tape' or 'Resin' or 'Undulation':
                for lay, d in enumerate(grids[self.layer:]):
                    if d[i].objectType == 'Polygon':
                        d[i].objectType = 'Tape'
                        d[i].layer = self.layer+lay+1
                        d[i].angle = self.angle
                        self.tapePath[index] = (i,d[i])
                        break
                    else: continue 
            else: continue

        # define undulation regions
        objByLayer = [[(ii,obj2) for (ii,obj2) 
                        in self.tapePath if obj2.layer == s] for s 
                        in range(1,len(grids)+1)]
        for n in range(len(objByLayer)-1):
            mergedLayer1 = cascaded_union(
                [obj5 for (obj5ID, obj5) in objByLayer[n]])
            layer2 = objByLayer[n+1]
            touches = [(objID, obj4) for (objID, obj4) 
                            in layer2 
                            if obj4.touches(mergedLayer1)]
            for (objid, obj5) in touches:
                    cellTop = grids[n+1][objid]
                    cellBottom = grids[n][objid]
                    if cellTop.objectType == 'Undulation':
                        cellTop.angle = self.angle
                    else:
                        cellTop.objectType = 'Undulation'
                        cellTop.angle = self.angle
                        # cellTop.layer = n
                    if cellBottom.objectType == 'Undulation':
                        cellBottom.angle = self.angle
                        self.tapePath.append([objid, cellBottom])
                    else:
                        cellBottom.objectType = 'Undulation'
                        cellBottom.angle = self.angle
                        cellBottom.layer = n+1
                        self.tapePath.append([objid, cellBottom])

        # create bordering resin regions
        self.tapeCentroid= bounds.centroid
        if self.angle == 90:
            # specified clockwise
            xx, yy = rotatedBounds.exterior.xy
            x0 = map(lambda m: round(m,6), xx)
            y0 = map(lambda n: round(n,6), yy)
            rCoords = zip(x0,y0)
            # left side
            resinTopRotated = Polygon([(rCoords[1][0]-resinW, rCoords[1][1]),
                        (rCoords[1][0], rCoords[1][1]),
                        (rCoords[0][0],rCoords[0][1]),
                        (rCoords[0][0]-resinW, rCoords[0][1])])
            # right side
            resinBottomRotated = Polygon([(rCoords[2][0], rCoords[2][1]),
                        (rCoords[2][0]+resinW, rCoords[2][1]),
                        (rCoords[3][0]+resinW,rCoords[3][1]),
                        (rCoords[3][0], rCoords[3][1])])
        else:
            resinTop = Polygon([(self.coords[0][0], self.coords[0][1]+resinW),
                        (self.coords[1][0], self.coords[1][1]+resinW),
                        (self.coords[1][0],self.coords[1][1]),
                        (self.coords[0][0], self.coords[0][1])])
            resinBottom = Polygon([(self.coords[2][0], self.coords[2][1]),
                        (self.coords[3][0], self.coords[3][1]),
                        (self.coords[3][0],self.coords[3][1]-resinW),
                        (self.coords[2][0], self.coords[2][1]-resinW)])
            resinTopRotated = affinity.rotate(
                resinTop, self.angle, self.tapeCentroid)
            resinBottomRotated = affinity.rotate(
                resinBottom, self.angle, self.tapeCentroid)
        resinBounds = cascaded_union(
            [resinTopRotated,resinBottomRotated]).buffer(
                1*10**-8, join_style=2)
        self.resinRegions = [(i,plygn) for (i,plygn) 
                                in grids[self.layer-1].iteritems() 
                                if plygn.within(resinBounds)]

        # lift up overlapping resin regions
        for index3 in range(len(self.resinRegions)):
            iii, obj3 = self.resinRegions[index3]
            if obj3.objectType == 'Polygon':
                obj3.objectType = 'Resin'
                obj3.angle = self.angle
            elif obj3.objectType == 'Tape' or 'Resin' or 'Undulation':
                for ly, dd in enumerate(grids[self.layer:]):
                    if dd[iii].objectType == 'Polygon':
                        dd[iii].objectType = 'Resin'
                        dd[iii].layer = self.layer+ly+1
                        dd[iii].angle = self.angle
                        break
                    else: continue 
    
    def __str__(self):
        return 'MachinePass: Angle = {}'.format(self.angle)

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    grids = createGrids(tapeAngles=(0,90), tapeWidths=(10,10))

    one = machinePass(grids).tapePath
    two = machinePass(grids, angle=90).tapePath

    print [x.layer for (unuse,x) in two if x.objectType == 'Undulation' ]
    objPlot(grids, (0,90), 'Tape')
    objPlot(grids, (0,90), 'Resin')
    objPlot(grids, (0,90), 'Undulation')
# ------------------------------------------------------------------------

    


    
