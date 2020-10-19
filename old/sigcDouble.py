from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from shapely.ops import cascaded_union, linemerge, unary_union, polygonize
from copy import deepcopy
from math import cos, radians, tan
from itertools import combinations

'''
This script is used to define the geometry of interlaced laminates.
This module is imported by the tapePlacement script, which is itself called
by the shapelyToAbaqus script.

Geometries are created using the Python Shapely library. The script first
creates multiple polygon grids, each representing a layer in a laminate.

The afpPass class can be used to place tapes to form a laminate.
'''

# ----------------------------------------------------------------------------


def createGrids(tapeAngles, tapeWidths, undulationWidth=1.0, sample=None):
    '''
    Function to partition the specimen into a grid of polygons depending on the
    tape angles of the laminate.

    The specimen region (150mm x 100mm) is separated into Shapely Polygons.
    The tapes are bounded by undulation/resin regions (half the width of the
    uw parameter)

    '''
    # define specimen boundaries
    if sample is None:
        sample = Polygon(
            [(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])

    xMax = sample.exterior.xy[0][2]
    yMax = sample.exterior.xy[1][2]

    # define mirror point (used to mirror the Polygon boundaries)
    mirrorPoint = Point([(0.0, 0.0), (0.0, 0.0)])

    nLayers = len(tapeAngles)  # number of layers in the laminate
    tapes = zip(tapeAngles, tapeWidths)
    partitionLines = []
    uw = undulationWidth
    for tape in tapes:
        a = tape[0]  # angle
        w = tape[1]  # width
        if a == 90:
            offset = w
            maxOffset = xMax - w / 2.0
            numberOffsets = int(maxOffset / offset)
            offsetList = [offset * k for k
                          in range(0, numberOffsets + 1)]
            for os in offsetList:
                tapeLineCoords = [((w / 2.0) - uw + os, 300.0),
                                  ((w / 2.0) - uw + os, -300.0)]
                resinLineCoords1 = ([((w / 2.0) + os, 300.0),
                                     ((w / 2.0) + os, -300.0)])
                resinLineCoords2 = ([((w / 2.0) + os + uw, 300.0),
                                     ((w / 2.0) + os + uw, -300.0)])
                tapeLine = LineString(tapeLineCoords)
                resinLine1 = LineString(resinLineCoords1)
                resinLine2 = LineString(resinLineCoords2)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])
            reflectedLines = [scale(line, xfact=-1, origin=mirrorPoint)
                              for line in partitionLines]
            partitionLines = partitionLines + reflectedLines

        else:
            offset = w / cos(radians(a))
            maxOffset = (yMax - (w / 2.0) * cos(radians(a))
                         + xMax * tan(radians(a)))
            numberOffsets = int(maxOffset / offset)
            offsetList = [offset * k for k in range(0, numberOffsets + 20)]
            for os in offsetList:
                tapeLineCoords = [(-300.0, (w / 2.0) - uw + os),
                                  (300.0, (w / 2.0) - uw + os)]
                resinLineCoords1 = ([(-300.0, (w / 2.0) + os),
                                     (300.0, (w / 2.0) + os)])
                resinLineCoords2 = ([(-300.0, (w / 2.0) + os + uw),
                                     (300.0, (w / 2.0) + os + uw)])
                rotPoint = Point([(0.0, os), (0.0, os)])
                tapeLine = rotate(LineString(tapeLineCoords), a, rotPoint)
                resinLine1 = rotate(LineString(resinLineCoords1), a, rotPoint)
                resinLine2 = rotate(LineString(resinLineCoords2), a, rotPoint)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])
            reflectedLines = [rotate(line, 180.0, origin=mirrorPoint)
                              for line in partitionLines]
            partitionLines = partitionLines + reflectedLines

    # collection of individual linestrings for splitting in a list and add
    # the polygon lines to it.
    partitionLines.append(sample.boundary)
    merged_lines = linemerge(partitionLines)
    border_lines = unary_union(merged_lines)
    decomposition = polygonize(border_lines)
    trimmed = [plygn for plygn in decomposition if plygn.within(sample)]
    polyDict = {key: value for (key, value) in enumerate(trimmed)}
    layerGrids = {l: deepcopy(polyDict) for l in range(1, nLayers + 1)}
    return layerGrids


# initialize additional polygon attributes
setattr(Polygon, 'objectType', None)
setattr(Polygon, 'angle', None)  # create attribute in Polygon class
setattr(Polygon, 'layer', 1)  # create attribute in Polygon class
setattr(Polygon, 'objectID', None)


class afpPass():
    '''
    Class defining a single pass of a an AFP machine
    inputs: grid of the specimen, angle of the pass, coords, width of resin
    NOTE: a pass should be specified using its coordinates when horizontal
    and then rotated. Do not rotate outside of the afpPass class
    '''

    def __init__(self, grids, angle=0, coords=None, undulationWidth=1.0):
        self.angle = angle  # assign to instance variable
        self.uw = undulationWidth  # assigned to uw for brevity, default value = 1
        # set tape dimensions or use default
        if coords is None:
            w = 30.0
            self.coords = ([(-300.0, (w / 2) - self.uw), (300.0, (w / 2) - self.uw),
                            (300.0, (-w / 2) + self.uw), (-300.0, (-w / 2) + self.uw)])
        else:
            self.coords = coords
        # define tape boundaries
        self.bounds = Polygon(self.coords)
        self.rotatedBounds = rotate(self.bounds,
                               self.angle).buffer(1 * 10**-8, join_style=2)
        # identify which Polygons from the createGrids function are within
        # the boundaries of the rotated tape
        self.tapePath = [(obj1ID, obj1) for (obj1ID, obj1)
                         in grids[1].iteritems()
                         if obj1.within(self.rotatedBounds)]

        self.raiseTapes(grids)
        self.createResinRegions(grids, self.uw)
        self.raiseResin(grids)
        self.createUndulations(grids)
        self.defineConnectivity()

    def createResinRegions(self, grids, uw):
        # create bordering resin regions
        self.tapeCentroid = self.bounds.centroid
        if self.angle == 90:

            # specified clockwise (1,2,3,4)
            # 1 _______________________________________2
            #  |______________________________________|   Resin region (top)
            # 4|                                      |3
            # 1|______________________________________|2  Tape region
            #  |______________________________________|   Resin region (bottom)
            # 3                                        4

            xx, yy = self.rotatedBounds.exterior.xy
            x0 = map(lambda m: round(m, 6), xx)
            y0 = map(lambda n: round(n, 6), yy)
            rCoords = zip(x0, y0)
            # left side
            resinTopRotated = Polygon([(rCoords[1][0] - uw, rCoords[1][1]),
                                       (rCoords[1][0], rCoords[1][1]),
                                       (rCoords[0][0], rCoords[0][1]),
                                       (rCoords[0][0] - uw, rCoords[0][1])])
            # right side
            resinBottomRotated = Polygon([(rCoords[2][0], rCoords[2][1]),
                                          (rCoords[2][0] + uw, rCoords[2][1]),
                                          (rCoords[3][0] + uw, rCoords[3][1]),
                                          (rCoords[3][0], rCoords[3][1])])
        else:
            resinTop = Polygon([(self.coords[0][0], self.coords[0][1] + uw),
                                (self.coords[1][0], self.coords[1][1] + uw),
                                (self.coords[1][0], self.coords[1][1]),
                                (self.coords[0][0], self.coords[0][1])])
            resinBottom = Polygon([(self.coords[2][0], self.coords[2][1]),
                                   (self.coords[3][0], self.coords[3][1]),
                                   (self.coords[3][0], self.coords[3][1] - uw),
                                   (self.coords[2][0], self.coords[2][1] - uw)])
            resinTopRotated = rotate(resinTop, self.angle, self.tapeCentroid)
            resinBottomRotated = rotate(resinBottom, self.angle,
                                        self.tapeCentroid)
        resinBounds = cascaded_union([resinTopRotated,
                                      resinBottomRotated]).buffer(1 * 10**-8,
                                                                  join_style=2)
        self.resinRegions = [(obj5ID, obj5) for (obj5ID, obj5)
                             in grids[1].iteritems()
                             if obj5.within(resinBounds)]

    def raiseResin(self, grids):
        # lift up overlapping resin regions
        for j in range(len(self.resinRegions)):
            obj4ID, obj4 = self.resinRegions[j]
            if obj4.objectType is None:
                setattr(obj4, 'objectType', 'Resin')
                self.setAngle(obj4, self.angle)
            elif obj4.objectType == 'Tape':
                k = 2
                while k <= len(grids):
                    objGrid2 = grids[k]  # next grid layer
                    if objGrid2[obj4ID].objectType is None:
                        setattr(objGrid2[obj4ID], 'objectType', 'Resin')
                        objGrid2[obj4ID].layer = k
                        self.setAngle(objGrid2[obj4ID], self.angle)
                        self.resinRegions[j] = (obj4ID, objGrid2[obj4ID])
                        break
                    # elif obj4.objectType == 'Undulation' and len(obj4.angle) == 1:
                    #     self.setAngle(obj4, self.angle)
                    #     print 'Triggered'
                    else:
                        k += 1

        self.tapePath.extend(self.resinRegions)
        self.inbetween(grids, 'Resin')

    def raiseTapes(self, grids):
        '''
        The following loop iterates over each Polygon identified as being
        within the bounds of the current machine pass.
        Each Polygon has an attribute called 'objectType'. The loop checks
        the Polygon's attribute type. If it is 'Polygon' that means that the
        object has not yet been identified as a Tape or Undulation region.
        In this case, the loop sets the objectType attribute of the Polygon
        to 'Tape'.
        If the objectType of a Polygon is 'Tape' or 'Undulation' that means
        that a previous pass has already assigned an objectType to the
        Polygon in question. In other words, there is already a tape placed
        in this area in this layer.
        If this is the case, the loop sets the attribute of the same Polygon
        in the next grid layer to 'Tape'. Essentially, if a tape crossing is
        detected then the Tape is placed in the next layer.
        '''
        for i in range(len(self.tapePath)):
            obj2ID, obj2 = self.tapePath[i]
            if obj2.objectType is None:
                setattr(obj2, 'objectType', 'Tape')
                self.setAngle(obj2, self.angle)
            elif obj2.objectType == 'Tape':
                # check layer above
                m = 2
                while m <= len(grids):
                    objGrid = grids[m]  # next grid layer
                    if objGrid[obj2ID].objectType is None:
                        setattr(objGrid[obj2ID], 'objectType', 'Tape')
                        objGrid[obj2ID].layer = m
                        self.setAngle(objGrid[obj2ID], self.angle)
                        self.tapePath[i] = (obj2ID, objGrid[obj2ID])
                        break
                    else:
                        m += 1
        self.inbetween(grids, 'Tape')

    def inbetween(self, grids, objType):
        objByLayer = [[(objID, obj) for (objID, obj)
                       in self.tapePath if obj.layer == s and obj.objectType == objType] for s
                      in range(1, len(grids) + 1)]
        for n in range(1, len(objByLayer)):
            # raise regions between objects in the same layer
            allObjInLayer = cascaded_union([x for (xID, x) in objByLayer[n]])
            try:
                allObjBuff = [obj.buffer(2 * self.uw, join_style=2)
                            for obj in allObjInLayer.geoms]
                allObjConvexDiff = cascaded_union(
                    [a.intersection(b) for a, b in combinations(allObjBuff, 2)])
                for (zID, z) in self.tapePath:
                    if z.within(allObjConvexDiff):
                        cell = grids[n + 1][zID]
                        cell.objectType = objType
                        cell.layer = n + 1
                        self.setAngle(cell, self.angle)
                        idx = [i for (i, el) in self.tapePath].index(zID)
                        self.tapePath[idx] = (zID, cell)
            except AttributeError:
                pass

    def createUndulations(self, grids):
        '''
        Create undulation regions
        '''

        objByLayer = [[(objID, obj) for (objID, obj)
                       in self.tapePath if obj.layer == s] for s
                      in range(1, len(grids) + 1)]
        for n in range(1, len(objByLayer)):
            topLayer = cascaded_union([x for (xID, x) in objByLayer[n]])
            bottomLayer = objByLayer[n - 1]
            topLayerBuff = topLayer.buffer(2 * self.uw + 1.0, join_style=2)
            bottomObjTouchingTopObj = [(zID, z) for (zID, z) in bottomLayer
                                       if z.within(topLayerBuff)]
            for (cellID, cell) in bottomObjTouchingTopObj:
                cellBottom = grids[n][cellID]
                cellTop = grids[n + 1][cellID]
                cellTop.objectType = 'Undulation'
                cellBottom.objectType = 'Undulation'
                if len(cellBottom.angle) != 2:
                    self.setAngle(cellBottom, self.angle)
                self.setAngle(cellTop, self.angle)
                self.tapePath.append([cellID, cellTop])

    def defineConnectivity(self):
        self.connectivity = []
        for (obj9ID, obj9) in self.tapePath:
            ident = '{}-{}'.format(obj9.layer, obj9ID)
            setattr(obj9, 'objectID', ident)
            self.connectivity.append(ident)

    def setAngle(self, obj, angle):
        '''
        Cannot setattr(Polygon, 'angle, set()) because sets are  mutable.
        This would result in all the Polygons sharing the same angle attribute.
        Instead initialize the angle attribute with None and then use setAngle
        to either create a set with the first angle or add to the set if an
        angle has already been assigned.
        '''
        if obj.angle is None:
            obj.angle = [angle, ]
        else:
            obj.angle.append(angle)

    def __str__(self):
        return 'afpPass: Angle = {}'.format(self.angle)

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.


if __name__ == '__main__':
    from objectPlot import objPlot
    import matplotlib.pyplot as plt

    tapeAng = (0, 90, 45)
    tapeW = (30.0, ) * len(tapeAng)

    # create grids
    grids = createGrids(tapeAngles=tapeAng, tapeWidths=tapeW,
                        undulationWidth=1.0)

    # place tapes
    # for ang in [0, 90, 45]:
        # passs = afpPass(grids, angle=ang, undulationWidth=1.0)

    passs = afpPass(grids, angle=0, undulationWidth=1.0)
    # objPlot(grids, tapeAng, 'Tape')
    # objPlot(grids, tapeAng, 'Resin')
    # objPlot(grids, tapeAng, 'Undulation')
    w = 30.0
    uw = 1.0
    crds = ([(-300.0, (w / 2) - uw -w), (300.0, (w / 2) - uw - w),
                    (300.0, (-w / 2) + uw - w), (-300.0, (-w / 2) + uw - w)])
    passs = afpPass(grids, coords=crds, angle=0, undulationWidth=1.0)
    passs = afpPass(grids, angle=90, undulationWidth=1.0)
    # objPlot(grids, tapeAng, 'Tape')
    # objPlot(grids, tapeAng, 'Resin')
    # objPlot(grids, tapeAng, 'Undulation')

    # passs = afpPass(grids, angle=45, undulationWidth=1.0)
    objPlot(grids, tapeAng, 'Tape')
    objPlot(grids, tapeAng, 'Resin')
    objPlot(grids, tapeAng, 'Undulation')
