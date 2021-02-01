from shapely.geometry import Polygon, Point, LineString, JOIN_STYLE
from shapely import affinity
from shapely.ops import cascaded_union, split, linemerge, unary_union, polygonize
from copy import deepcopy
import numpy as np
import math
from objectPlot import objPlot

'''
This script is used to define the geometry of interlaced laminates.
This module is imported by the tapePlacement script, which is itself called
by the shapelyToAbaqus script. 
Geometries are created using the Python Shapely library. The script first
creates multiple polygon grids, each representing a layer in a laminate. 
A machinepass is defined by the width of a tape and the 
'''

# ----------------------------------------------------------------------------

# Function to partition the specimen into a grid of polygons depending on the 
# tape angles of the laminate
def createGrids(tapeAngles, tapeWidths, undulationWidth=1.0, sample=None):
    # define specimen boundaries
    if sample == None:
        sample = Polygon(
            [(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])
    
    nLayers = len(tapeAngles) # number of layers in the laminate 
    tapes = zip(tapeAngles, tapeWidths)
    partitionLines = []
    uw = undulationWidth
    for tape in tapes:
        a = tape[0]  # angle
        w = tape[1]  # width
        if a == 90:
            offset = w
            maxOffset = 50.0-w/2.0
            numberOffsets = int(maxOffset/offset)
            offsetList = [offset*k for k 
                            in range(-numberOffsets-1,numberOffsets+1)]
            for os in offsetList:
                tapeLineCoords = [((w/2.0)-uw+os, 100.0), ((w/2.0)-uw+os, -100.0)]
                resinLineCoords1 = (
                        [((w/2.0)+os, 100.0),
                        ((w/2.0)+os, -100.0)])
                resinLineCoords2 = (
                        [((w/2.0)+os+uw, 100.0),
                        ((w/2.0)+os+uw, -100.0)])

                tapeLine = LineString(tapeLineCoords)
                resinLine1 = LineString(resinLineCoords1)
                resinLine2 = LineString(resinLineCoords2)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])

        else:
            offset = w/math.cos(math.radians(a))
            maxOffset = (75.0 - (w/2.0)*math.cos(math.radians(a))
                    + 50.0*math.tan(math.radians(a)))
            numberOffsets = int(maxOffset/offset)
            offsetList = [offset*k for k 
                            in range(-numberOffsets-5,numberOffsets+5)]
            for os in offsetList:
                tapeLineCoords = [(-100.0, (w/2.0)-uw+os),
                                  (100.0, (w/2.0)-uw+os)]
                resinLineCoords1 = (
                        [(-100.0, (w/2.0)+os), (100.0, (w/2.0)+os)])
                resinLineCoords2 = (
                        [(-100.0, (w/2.0)+os+uw), (100.0, (w/2.0)+os+uw)])
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
    layerGrids = {l: deepcopy(polyDict) for l in range(1,nLayers+1)}
    return layerGrids

# initialize additional polygon attributes
setattr(Polygon, 'objectType', None)
setattr(Polygon, 'angle', None) # create attribute in Polygon class
setattr(Polygon, 'layer', 1) # create attribute in Polygon class
setattr(Polygon, 'objectID', None)

# Class defining a single pass of a an ATP machine
# inputs: grid of the specimen, angle of the pass, coords, width of resin
# NOTE: a pass should be specified using its coordinates when horizontal 
# and then rotated. Do not rotate outside of the machinePass class
class machinePass():
    def __init__(self, grids, angle=0, coords=None, undulationWidth=1.0):
        self.angle = angle
        uw = undulationWidth
        if coords == None:
            w = 10.0
            self.coords = (
                    [(-100.0, (w/2)-uw), (100.0, (w/2)-uw),
                     (100.0,(-w/2)+uw), (-100.0, (-w/2)+uw)])
        else:
            self.coords = coords
        
        # define tape boundaries
        bounds = Polygon(self.coords)
        rotatedBounds = affinity.rotate(
                bounds, self.angle).buffer(1*10**-8, join_style=2)

        # identify which Polygons from the createGrids function are within
        # the boundaries of the rotated tape
        self.tapePath = [(obj1ID,obj1) for (obj1ID,obj1) 
                            in grids[1].iteritems() 
                            if obj1.within(rotatedBounds)]

        # The following loop iterates over each Polygon identified as being
        # within the bounds of the current machine pass. 
        # Each Polygon has an attribute called 'objectType'. The loop checks
        # the Polygon's attribute type. If it is 'Polygon' that means that the 
        # object has not yet been identified as a Tape or Undulation region. 
        # In this case, the loop sets the objectType attribute of the Polygon
        # to 'Tape'.
        # If the objectType of a Polygon is 'Tape' or 'Undulation' that means
        # that a previous pass has already assigned an objectType to the 
        # Polygon in question. In other words, there is already a tape placed
        # in this area in this layer. 
        # If this is the case, the loop sets the attribute of the same Polygon
        # in the next grid layer to 'Tape'. Essentially, if a tape crossing is
        # detected then the Tape is placed in the next layer.

        for i in range(len(self.tapePath)):
            obj2ID, obj2 = self.tapePath[i]
            if obj2.objectType == None:
                setattr(obj2,'objectType','Tape')
                self.setAngle(obj2,self.angle)
            elif obj2.objectType == 'Tape' or 'Undulation' or 'Resin':
                # check layer above
                m = 2
                while m <= len(grids):                    
                    objGrid = grids[m] # next grid layer
                    if objGrid[obj2ID].objectType == None:
                        setattr(objGrid[obj2ID],'objectType','Tape')
                        objGrid[obj2ID].layer = m
                        self.setAngle(objGrid[obj2ID],self.angle)
                        self.tapePath[i] = (obj2ID,objGrid[obj2ID])
                        break
                    else:
                        m += 1

        # create bordering resin regions
        self.tapeCentroid= bounds.centroid
        if self.angle == 90:
            # specified clockwise
            xx, yy = rotatedBounds.exterior.xy
            x0 = map(lambda m: round(m,6), xx)
            y0 = map(lambda n: round(n,6), yy)
            rCoords = zip(x0,y0)
            # left side
            resinTopRotated = Polygon([(rCoords[1][0]-uw, rCoords[1][1]),
                        (rCoords[1][0], rCoords[1][1]),
                        (rCoords[0][0],rCoords[0][1]),
                        (rCoords[0][0]-uw, rCoords[0][1])])
            # right side
            resinBottomRotated = Polygon([(rCoords[2][0], rCoords[2][1]),
                        (rCoords[2][0]+uw, rCoords[2][1]),
                        (rCoords[3][0]+uw,rCoords[3][1]),
                        (rCoords[3][0], rCoords[3][1])])
        else:
            resinTop = Polygon([(self.coords[0][0], self.coords[0][1]+uw),
                        (self.coords[1][0], self.coords[1][1]+uw),
                        (self.coords[1][0],self.coords[1][1]),
                        (self.coords[0][0], self.coords[0][1])])
            resinBottom = Polygon([(self.coords[2][0], self.coords[2][1]),
                        (self.coords[3][0], self.coords[3][1]),
                        (self.coords[3][0],self.coords[3][1]-uw),
                        (self.coords[2][0], self.coords[2][1]-uw)])
            resinTopRotated = affinity.rotate(
                resinTop, self.angle, self.tapeCentroid)
            resinBottomRotated = affinity.rotate(
                resinBottom, self.angle, self.tapeCentroid)
        resinBounds = cascaded_union(
            [resinTopRotated,resinBottomRotated]).buffer(
                1*10**-8, join_style=2)
        self.resinRegions = [(obj7ID,obj7) for (obj7ID,obj7) 
                                in grids[1].iteritems() 
                                if obj7.within(resinBounds)]
        

        # lift up overlapping resin regions
        for j in range(len(self.resinRegions)):
            obj8ID, obj8 = self.resinRegions[j]
            if obj8.objectType == None:
                setattr(obj8,'objectType','Resin')
                self.setAngle(obj8,self.angle)
            elif obj8.objectType == 'Tape' or 'Resin' or 'Undulation':
                k = 2 
                while k <= len(grids):
                    objGrid2 = grids[k] # next grid layer
                    if objGrid2[obj8ID].objectType == None:
                        setattr(objGrid2[obj8ID],'objectType','Resin')
                        objGrid2[obj8ID].layer = k
                        self.setAngle(objGrid2[obj8ID],self.angle)
                        self.resinRegions[j] = (obj8ID,objGrid2[obj8ID])
                        break
                    else: 
                        k += 1

        self.tapePath.extend(self.resinRegions)

        # define undulation regions
        objByLayer = [[(obj3ID,obj3) for (obj3ID,obj3) 
                        in self.tapePath if obj3.layer == s] for s 
                        in range(1,len(grids)+1)]
        for n in range(1, len(objByLayer)):
            bottomLayer = cascaded_union(
                [obj4 for (obj4ID, obj4) in objByLayer[n-1]])
            topLayer = objByLayer[n]
            topObjTouchingBottomObj = [(obj5ID, obj5) for (obj5ID, obj5) 
                            in topLayer 
                            if obj5.touches(bottomLayer)]
            for (obj6ID, obj6) in topObjTouchingBottomObj:
                    cellTop = grids[n+1][obj6ID]
                    cellTop.objectType = 'Undulation'
                    self.setAngle(cellTop,self.angle)
                    cellBottom = grids[n][obj6ID]
                    cellBottom.objectType = 'Undulation'
                    self.setAngle(cellBottom,self.angle)
                    self.tapePath.append([obj6ID, cellBottom])

        self.connectivity = []
        for (objID,obj) in self.tapePath:
            ident = '{}-{}'.format(obj.layer,objID) 
            setattr(obj,'objectID', ident)
            self.connectivity.append(ident)

        objPlot(grids, (0,90,45), 'Undulation')
        objPlot(grids, (0,90,45), 'Resin')

        # print self.connectivity
        # print len(self.connectivity)
    # Cannot setattr(Polygon, 'angle, set()) because sets are  mutable.
    # This would result in all the Polygons sharing the same angle attribute.
    # Instead initialize the angle attribute with None and then use setAngle
    # to either create a set with the first angle or add to the set if an 
    # angle has already been assigned. 
    def setAngle(self, obj, angle): 
        if obj.angle == None:
            # obj.angle = set((angle, ))
            obj.angle = angle
        else:
            # obj.angle = set((angle, ))
            obj.angle = angle
            # obj.angle.add(angle)

    def __str__(self):
        return 'MachinePass: Angle = {}'.format(self.angle)

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (obj1ID.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from objectPlot import objPlot

    tapeAng = (0,90,45)
    grids = createGrids(tapeAngles=tapeAng, tapeWidths=(10,10,10),
            undulationWidth=1.0)

    for (key,poly) in grids[1].iteritems():
        x,y = poly.exterior.xy
        plt.plot(x,y)
    plt.grid(True)
    plt.show()

    for ang in tapeAng:
        passs = machinePass(grids, angle=ang, undulationWidth=1.0)

    # for key, polyObj in grids[1].iteritems():
    #     if polyObj.objectType == 'Undulation':
    #         print polyObj.angle

    objPlot(grids, tapeAng, 'Tape')
    objPlot(grids, tapeAng, 'Resin')
    objPlot(grids, tapeAng, 'Undulation')
# ------------------------------------------------------------------------
