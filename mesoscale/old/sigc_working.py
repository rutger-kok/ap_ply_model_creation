from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from shapely.ops import cascaded_union, linemerge, unary_union, polygonize
from copy import deepcopy
from math import cos, radians, tan

'''
This script is used to define the geometry of interlaced laminates.
This module is imported by the tapePlacement script, which is itself called
by the shapelyToAbaqus script.

Geometries are created using the Python Shapely library. The script first
creates multiple polygon grids, each representing a layer in a laminate.

The machinePass class can be used to place tapes to form a laminate.
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


# Class defining a single pass of a an ATP machine
# inputs: grid of the specimen, angle of the pass, coords, width of resin
# NOTE: a pass should be specified using its coordinates when horizontal
# and then rotated. Do not rotate outside of the machinePass class
class machinePass():
    def __init__(self, grids, angle=0, coords=None, undulationWidth=1.0):

        self.angle = angle  # assign to instance variable
        uw = undulationWidth  # assigned to uw for brevity, default value = 1

        # default tape dimensions (used unless other coords are supplied)
        if coords is None:
            w = 10.0
            self.coords = ([(-300.0, (w / 2) - uw), (300.0, (w / 2) - uw),
                            (300.0, (-w / 2) + uw), (-300.0, (-w / 2) + uw)])
        else:
            self.coords = coords

        # define tape boundaries
        bounds = Polygon(self.coords)
        rotatedBounds = rotate(bounds,
                               self.angle).buffer(1 * 10**-8, join_style=2)

        # identify which Polygons from the createGrids function are within
        # the boundaries of the rotated tape
        self.tapePath = [(obj1ID, obj1) for (obj1ID, obj1)
                         in grids[1].iteritems()
                         if obj1.within(rotatedBounds)]

        # The following loop iterates over each Polygon identified as being
        # within the bounds of the current machine pass.
        # Each Polygon has an attribute called 'objectType'. The loop checks
        # the Polygon's attribute type. If it is 'None' that means that the
        # object has not yet been identified as a Tape or Undulation region.
        # In this case, the loop sets the objectType attribute of the Polygon
        # to 'Tape'.
        # If the objectType of a Polygon is 'Tape' that means
        # that a previous pass has already assigned an objectType to the
        # Polygon in question. In other words, there is already a tape placed
        # in this area in this layer.
        # If this is the case, the loop sets the attribute of the same Polygon
        # in the next grid layer to 'Tape'. Essentially, if a tape crossing is
        # detected then the Tape is placed in the next layer.

        for i in range(len(self.tapePath)):
            obj2ID, obj2 = self.tapePath[i]
            if obj2.objectType is None:
                setattr(obj2, 'objectType', 'Tape')
                self.setAngle(obj2, self.angle)
            elif obj2.objectType == 'Tape' or 'Resin' or 'Undulation':
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

        # define undulation regions
        objByLayer = [[(obj5ID, obj5) for (obj5ID, obj5)
                      in self.tapePath if obj5.layer == s] for s
                      in range(1, len(grids) + 1)]
        for n in range(1, len(objByLayer)):
            bottomLayer = cascaded_union(
                [obj6 for (obj6ID, obj6) in objByLayer[n - 1]])
            bottomLayerBuff = bottomLayer.buffer(undulationWidth + 0.0001,
                                                 join_style=2)
            topLayer = objByLayer[n]
            topObjTouchingBottomObj = [(obj7ID, obj7) for (obj7ID, obj7)
                                       in topLayer
                                       if obj7.within(bottomLayerBuff)]
            for (obj8ID, obj8) in topObjTouchingBottomObj:
                    cellBottom = grids[n][obj8ID]
                    cellTop = grids[n + 1][obj8ID]
                    cellTopAngle = [cellBottom.angle[-1]]
                    cellTop.objectType = 'Undulation'
                    cellBottom.objectType = 'Undulation'
                    if len(cellBottom.angle) == 2:
                        cellBottom.angle = [cellBottom.angle[-1], self.angle]
                    else:
                        self.setAngle(cellBottom, self.angle)
                    self.setAngle(cellTop, cellTopAngle)
                    self.tapePath.append([obj8ID, cellBottom])

        self.connectivity = []
        for (obj9ID, obj9) in self.tapePath:
            ident = '{}-{}'.format(obj9.layer, obj9ID)
            setattr(obj9, 'objectID', ident)
            self.connectivity.append(ident)

        # create bordering resin regions
        self.tapeCentroid = bounds.centroid
        if self.angle == 90:

            # specified clockwise (1,2,3,4)
            # 1 _______________________________________2
            #  |______________________________________|   Resin region (top)
            # 4|                                      |3
            # 1|______________________________________|2  Tape region
            #  |______________________________________|   Resin region (bottom)
            # 3                                        4

            xx, yy = rotatedBounds.exterior.xy
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

        # lift up overlapping resin regions
        for j in range(len(self.resinRegions)):
            obj4ID, obj4 = self.resinRegions[j]
            if obj4.objectType is None:
                setattr(obj4, 'objectType', 'Resin')
                self.setAngle(obj4, self.angle)
            elif obj4.objectType == 'Tape' or 'Resin' or 'Undulation':
                k = 2
                while k <= len(grids):
                    objGrid2 = grids[k]  # next grid layer
                    if objGrid2[obj4ID].objectType is None:
                        setattr(objGrid2[obj4ID], 'objectType', 'Resin')
                        objGrid2[obj4ID].layer = k
                        self.setAngle(objGrid2[obj4ID], self.angle)
                        self.resinRegions[j] = (obj4ID, objGrid2[obj4ID])
                        break
                    else:
                        k += 1

        # self.tapePath.extend(self.resinRegions)

    # Cannot setattr(Polygon, 'angle, set()) because sets are  mutable.
    # This would result in all the Polygons sharing the same angle attribute.
    # Instead initialize the angle attribute with None and then use setAngle
    # to either create a set with the first angle or add to the set if an
    # angle has already been assigned.
    def setAngle(self, obj, angles):
        if obj.angle is None:
            if isinstance(angles, list):
                obj.angle = angles
            else:
                obj.angle = [angles]
        else:
            if isinstance(angles, list):
                obj.angle.extend(angles)
            else:
                obj.angle.append(angles)

    def __str__(self):
        return 'MachinePass: Angle = {}'.format(self.angle)

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.


if __name__ == '__main__':
    from objectPlot import objPlot

    tapeAng = (0, 90)
    tapeW = (25.0, ) * len(tapeAng)

    # create grids
    grids = createGrids(tapeAngles=tapeAng, tapeWidths=tapeW,
                        undulationWidth=1.0)

    # place tapes
    for ang in [0, 90]:
        passs = machinePass(grids, angle=ang, undulationWidth=1.0)

    # objPlot(grids, tapeAng, 'Tape')
    # objPlot(grids, tapeAng, 'Resin')
    objPlot(grids, tapeAng, 'Undulation')
