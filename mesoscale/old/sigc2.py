from shapely.geometry import Polygon, Point, LineString, shape, JOIN_STYLE, MultiLineString
from shapely import affinity
from shapely.ops import cascaded_union, split, linemerge, unary_union, polygonize
from shapely.strtree import STRtree
from itertools import count
from copy import deepcopy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

def str_to_class(classname):
    return getattr(sys.modules[__name__], classname)

sample = Polygon([(-50.0, -75.1), (50.0, -75.1), (50.0, 75.1), (-50.0, 75.1)])

w = 10.0
t = 0.2

setattr(Polygon, 'objectType', 'Polygon') # create attribute in Polygon class

def createGrids(tapeAngles, tapeWidths, r=1.0, sample=None):
    # define specimen boundaries
    if sample == None:
        sample = Polygon([(-50.0, -75.1), (50.0, -75.1), (50.0, 75.1), (-50.0, 75.1)])
    
    tapes = zip(tapeAngles, tapeWidths)
    partitionLines = []
    for tape in tapes:
        a = tape[0]  # angle
        w = tape[1]  # width
        if a == 90:
            offset = (w+2*r)
            maxOffset = 50.0-w/2.0
            numberOffsets = int(maxOffset/offset)
            offsetList = [offset*k for k in range(-numberOffsets-1,numberOffsets+1)]
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
            offsetList = [offset*k for k in range(-numberOffsets-1,numberOffsets+1)]
            for os in offsetList:
                tapeLineCoords = [(-100.0, (w/2.0)+os), (100.0, (w/2.0)+os)]
                resinLineCoords1 = (
                        [(-100.0, (w/2.0)+os+r), (100.0, (w/2.0)+os+r)])
                resinLineCoords2 = (
                        [(-100.0, (w/2.0)+os+2*r), (100.0, (w/2.0)+os+2*r)])
                rotPoint = Point([(0.0, os), (0.0, os)])
                tapeLine = affinity.rotate(LineString(tapeLineCoords), a, rotPoint)
                resinLine1 = affinity.rotate(LineString(resinLineCoords1), a, rotPoint)
                resinLine2 = affinity.rotate(LineString(resinLineCoords2), a, rotPoint)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])

    partitionLines.append(sample.boundary) # collection of individual linestrings for splitting in a list and add the polygon lines to it.
    merged_lines = linemerge(partitionLines)
    border_lines = unary_union(merged_lines)
    decomposition = polygonize(border_lines)     
    trimmed = [plygn for plygn in decomposition if plygn.within(sample)]
    layerGrids = [deepcopy(trimmed) for l in range(len(tapeAngles))]
    return layerGrids

grids = createGrids(tapeAngles=(0,90), tapeWidths=(10,10))

# for poly in grids[0]:
#     x,y = poly.exterior.xy
#     plt.plot(x,y)
#     plt.grid(True)
# plt.show()


class machinePass():
    def __init__(self, grids, layerNum=1, angle=0, coords=None, resinW=1.0):
        self.layer = layerNum
        self.angle = angle
        if coords == None:
            self.coords = (
                    [(-100.0, w/2), (100.0, w/2), (100.0, -w/2), (-100.0, -w/2)]
                    )
        else:
            self.coords = coords

        bounds = Polygon(self.coords)
        rotatedBounds = affinity.rotate(bounds, self.angle).buffer(1*10**-8, join_style=2)

        tapeRegions = [plygn for plygn in grids[self.layer-1] if plygn.within(rotatedBounds)]
        self.tapePath = sorted(
                    tapeRegions,
                    key = lambda x: (round(x.centroid.coords[0][0],6),
                        round(x.centroid.coords[0][1],6)))

        for index in range(len(self.tapePath)):
            obj = self.tapePath[index]
            x,y = obj.exterior.xy
            objCoords = zip(x,y)
            if obj.objectType == 'Polygon':
                newTape = Tape(coords=objCoords,angle=self.angle)
                grids[self.layer-1].append(newTape)
                grids[self.layer-1].remove(obj)
                self.tapePath[index] = newTape
            elif obj.objectType == 'Tape' or 'Resin':
                newObj = self.moveUp(obj, self.layer, 'Tape')
                self.tapePath[index] = newObj
            else: continue

        for index2 in range(len(self.tapePath)-1):
            obj2 = self.tapePath[index2]
            x2,y2 = obj2.exterior.xy
            obj2Coords = zip(x2,y2)
            neighbor1 = self.tapePath[index2-1]
            neighbor2 = self.tapePath[index2+1]
            if obj2.layer == self.layer and neighbor1.layer != neighbor2.layer:
                newUnd1 = Undulation(coords=obj2Coords,angle=self.angle, layer=neighbor1.layer)
                grids[self.layer-1].append(newUnd1)  # in same layer
                newUnd2 = Undulation(coords=obj2Coords,angle=self.angle, layer=neighbor2.layer+1)
                grids[self.layer].append(newUnd2)  # in next layer
                tempObj = [plygn for plygn in grids[self.layer] if plygn.within(obj2.buffer(1*10**-8))][0]
                grids[self.layer].remove(tempObj)  # in next layer
                grids[self.layer-1].remove(obj2)
                # self.tapePath[index2:index2+1] = (newUnd1, newUnd2)


        # create resin regions
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
            resinTopRotated = affinity.rotate(resinTop, self.angle, self.tapeCentroid)
            resinBottomRotated = affinity.rotate(resinBottom, self.angle, self.tapeCentroid)
        resinBounds = cascaded_union([resinTopRotated,resinBottomRotated]).buffer(1*10**-8, join_style=2)
        self.resinRegions = [plygn for plygn in grids[self.layer-1] if plygn.within(resinBounds)]

        for index3 in range(len(self.resinRegions)):
            obj3 = self.resinRegions[index3]
            x3,y3 = obj3.exterior.xy
            obj3Coords = zip(x3,y3)
            if obj3.objectType == 'Polygon':
                grids[self.layer-1].append(Resin(coords=obj3Coords,angle=self.angle))
                grids[self.layer-1].remove(obj3)
            elif obj3.objectType == 'Tape' or 'Resin' or 'Undulation':
                newObj = self.moveUp(obj3, self.layer, 'Resin')
                self.resinRegions[index3] = newObj
    
    def moveUp(self, obj, layerNumber, currentType):
        x,y = obj.exterior.xy
        objCoords = zip(x,y)
        a = True
        while a == True:
            n = 0
            cellAbove = [plygn for plygn in grids[layerNumber+n] if plygn.within(obj.buffer(1*10**-8, join_style=2))]
            if len(cellAbove)>1: print 'Error: multiple objects within'
            check = cellAbove[0]
            if check.objectType != 'Polygon':
                n += 1
            else:
                a = False
                class_ = str_to_class(currentType)
                newObj = class_(coords=objCoords, angle=self.angle, layer=layerNumber+1)
                grids[layerNumber+n].append(newObj)
                grids[layerNumber+n].remove(check)
        return newObj



# defines the Tape class, which is a subclass of the Shapely Polygon class.
class Tape(Polygon):
    _ids = count(0)  # counts the number of instances created
    _instances = []  # initializes a list of Tape instances 
    def __init__(self, coords=None, angle=0, layer=1):
        self.objectType = 'Tape'
        self.angle = angle
        self.cds = coords
        self.layer = layer
        Polygon.__init__(self,coords)

# the Resin class is a subclass of the Tape class. Aside from object creation, 
# the class defines its own __init__ method as its intersection behaviour is
# different (resin intersecting tape -> no undulation). 
class Resin(Polygon):
    _ids = count(0)
    _instances = []
    def __init__(self, coords=None, angle=0, layer=1):
        self.objectType = 'Resin'
        self.angle = angle
        self.cds = coords
        self.layer = layer
        Polygon.__init__(self,coords) 

# the Undulation class is a subclass of the Tape class. The Undulation class
# defines its own __init__ method which defines how the undulation regions
# are recorded.
class Undulation(Polygon):
    _ids = count(0)
    _instances = []
    def __init__(self, coords=None, angle=0, layer=1):
        self.objectType = 'Undulation'
        self.angle = angle
        self.cds = coords
        self.layer = layer
        Polygon.__init__(self,coords)

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.

if __name__ == '__main__':
    machinePass(grids)
    # machinePass(grids, coords=[(-100.0, -7.0), (100.0, -7.0), (100.0, -17.0), (-100.0, -17.0)])
    machinePass(grids, angle=90)
    # machinePass(grids, angle=45)

    print len(grids[0])
    print len(grids[1])
    # ------------------------------------------------------------------------

    # f1, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    # f1.suptitle('Undulation per Layer')

    # g = 0
    # for i in range(3):
    #     for j in range(3):
    #         if g <= len(grids)-1:
    #             tapesInLayer = [tape for tape in grids[g] if tape.objectType == 'Undulation']
    #             for t in tapesInLayer:
    #                 xi,yi = t.exterior.xy
    #                 axes[i,j].plot(xi,yi)
    #                 axes[i,j].set_title('Layer: {}'.format(g))
    #             xb, yb = sample.exterior.xy
    #             axes[i,j].plot(xb,yb)
    #             g += 1
    #         else:
    #             break

    # f2, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    # f2.suptitle('Tape per Layer')
    # g = 0
    # for i in range(3):
    #     for j in range(3):
    #         if g <= len(grids)-1:
    #             tapesInLayer = [tape for tape in grids[g] if tape.objectType == 'Tape']
    #             for t in tapesInLayer:
    #                 xi,yi = t.exterior.xy
    #                 axes[i,j].plot(xi,yi)
    #                 axes[i,j].set_title('Layer: {}'.format(g))
    #             xb, yb = sample.exterior.xy
    #             axes[i,j].plot(xb,yb)
    #             g += 1
    #         else:
    #             break

    # f3, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    # f3.suptitle('Resin per Layer')

    # g = 0
    # for i in range(3):
    #     for j in range(3):
    #         if g <= len(grids)-1:
    #             tapesInLayer = [tape for tape in grids[g] if tape.objectType == 'Resin']
    #             for t in tapesInLayer:
    #                 xi,yi = t.exterior.xy
    #                 axes[i,j].plot(xi,yi)
    #                 axes[i,j].set_title('Layer: {}'.format(g))
    #             xb, yb = sample.exterior.xy
    #             axes[i,j].plot(xb,yb)
    #             g += 1
    #         else:
    #             break

    # plt.show(f2)
    # plt.show(f3)
    
