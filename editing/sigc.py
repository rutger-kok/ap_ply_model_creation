from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from shapely.ops import cascaded_union, unary_union, polygonize
from copy import deepcopy
from math import cos, radians, tan
from itertools import combinations
from collections import OrderedDict

'''
This script is used to define the geometry of interlaced laminates.
This module is imported by the tapePlacement script, which is itself called
by the shapelyToAbaqus script.

Geometries are created using the Python Shapely library. The script first
creates multiple polygon grids, each representing a layer in a laminate.

The afpPass class can be used to place tapes to form a laminate.
'''

# ----------------------------------------------------------------------------

# initialize additional polygon attributes
setattr(Polygon, 'objectType', None)
setattr(Polygon, 'angle', None)  # create attribute in Polygon class
setattr(Polygon, 'layer', 1)  # create attribute in Polygon class
setattr(Polygon, 'objectID', None)


class Interlaced():
    def __init__(self, tapeAngles, tapeWidths, undulationWidth=1.0,
                 specimen=None):
        self.t_angles = tapeAngles
        self.t_widths = tapeWidths
        self.u_width = undulationWidth

        # define specimen boundaries
        if specimen is None:
            self.specimen = Polygon(
                [(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])
        else:
            self.specimen = specimen

        self.specimenBuffered = self.specimen.buffer(2.5, join_style=2)
        self._grids = self.createGrids()

    def makePass(self, passAngle=0, passCoords=None):
        '''
        Function defining a single pass of a an AFP machine
        inputs: grid of the specimen, angle of the pass, coords, width of resin
        NOTE: a pass should be specified using its coordinates when horizontal
        and then rotated. Do not rotate outside of the afpPass class
        '''

        # set tape dimensions or use default
        if passCoords is None:
            w = 30.0
            passCoords = ([(-300.0, (w / 2) - self.u_width),
                           (300.0, (w / 2) - self.u_width),
                           (300.0, (-w / 2) + self.u_width),
                           (-300.0, (-w / 2) + self.u_width)])

        # define tape boundaries
        passBounds = Polygon(passCoords)
        passRotBounds = rotate(
            passBounds, passAngle).buffer(1 * 10**-8, join_style=2)
        # identify which Polygons from the createGrids function are within
        # the boundaries of the rotated tape
        passPath = [(objID, obj) for (objID, obj)
                    in self._grids[1].iteritems() if obj.within(passRotBounds)]

        passPath = self.raiseTapes(passPath, passAngle)
        passPath = self.createUndulations(passPath, passAngle)
        self.createResinRegions(passPath, passAngle)
        self.grids = self.returnTrimmedGrid()
        connectivity = self.defineConnectivity(passPath)
        return connectivity

    def createGrids(self):
        '''
        Function to partition the specimen into a grid of polygons depending
        on the tape angles of the laminate.
        The specimen region (150mm x 100mm) is separated into Shapely Polygons.
        The tapes are bounded by undulation/resin regions (half the width of
        the u_width parameter)
        '''

        xMax = self.specimenBuffered.exterior.xy[0][2]
        yMax = self.specimenBuffered.exterior.xy[1][2]

        # define mirror point (used to mirror the Polygon boundaries)
        mirrorPoint = Point([(0.0, 0.0), (0.0, 0.0)])

        nLayers = len(self.t_angles)  # number of layers in the laminate
        tapes = zip(self.t_angles, self.t_widths)
        partitionLines = []
        for (a, w) in tapes:
            if a == 90:
                offset = w
                maxOffset = xMax - w / 2.0
                numOffsets = int(maxOffset / offset)
                offsetList = [offset * k for k in range(0, numOffsets + 1)]
                for ofs in offsetList:
                    tLineCoords = [((w / 2.0) - self.u_width + ofs, 300.0),
                                   ((w / 2.0) - self.u_width + ofs, -300.0)]
                    rLineCoords1 = ([((w / 2.0) + ofs, 300.0),
                                     ((w / 2.0) + ofs, -300.0)])
                    rLineCoords2 = ([((w / 2.0) + ofs + self.u_width, 300.0),
                                     ((w / 2.0) + ofs + self.u_width, -300.0)])
                    tLine = LineString(tLineCoords)
                    rLine1 = LineString(rLineCoords1)
                    rLine2 = LineString(rLineCoords2)
                    partitionLines.extend([tLine, rLine1, rLine2])
                reflectedLines = [scale(line, xfact=-1, origin=mirrorPoint)
                                  for line in partitionLines]
                partitionLines = partitionLines + reflectedLines

            else:
                offset = w / cos(radians(a))
                maxOffset = (yMax - (w / 2.0) * cos(radians(a))
                             + xMax * tan(radians(a)))
                numOffsets = int(maxOffset / offset)
                offsetList = [offset * k for k in range(0, numOffsets + 20)]
                for ofs in offsetList:
                    tLineCoords = [(-300.0, (w / 2.0) - self.u_width + ofs),
                                   (300.0, (w / 2.0) - self.u_width + ofs)]
                    rLineCoords1 = ([(-300.0, (w / 2.0) + ofs),
                                     (300.0, (w / 2.0) + ofs)])
                    rLineCoords2 = ([(-300.0, (w / 2.0) + ofs + self.u_width),
                                     (300.0, (w / 2.0) + ofs + self.u_width)])
                    rotPoint = Point([(0.0, ofs), (0.0, ofs)])
                    tLine = rotate(LineString(tLineCoords), a, rotPoint)
                    rLine1 = rotate(LineString(rLineCoords1), a, rotPoint)
                    rLine2 = rotate(LineString(rLineCoords2), a, rotPoint)
                    partitionLines.extend([tLine, rLine1, rLine2])
                reflectedLines = [rotate(line, 180.0, origin=mirrorPoint)
                                  for line in partitionLines]
                partitionLines = partitionLines + reflectedLines

        # collection of individual linestrings for splitting in a list and add
        # the polygon lines to it.
        partitionLines.append(self.specimenBuffered.boundary)
        partitionLines.append(self.specimen.boundary)
        border_lines = unary_union(partitionLines)
        decomposition = polygonize(border_lines)
        trimmed = [pgn for pgn in decomposition
                   if pgn.within(self.specimenBuffered)]
        polyDict = {key: value for (key, value) in enumerate(trimmed)}
        layerGrids = {l: deepcopy(polyDict) for l in range(1, nLayers + 1)}
        return layerGrids

    def returnTrimmedGrid(self):
        trimmedGrids = {}
        for layerNumber, gridLayer in self._grids.iteritems():
            trimmedGrids[layerNumber] = {}
            for (objID, obj) in gridLayer.iteritems():
                if obj.within(self.specimen):
                    trimmedGrids[layerNumber][objID] = obj
                else:
                    continue
        return trimmedGrids

    def createResinRegions(self, passPath, passAngle):
        objByLayer = [
            [(objID, obj) for (objID, obj) in passPath if obj.layer == s]
            for s in range(1, len(self._grids) + 1)]
        for n in range(1, len(objByLayer) + 1):
            allLayerObj = cascaded_union([x for (xID, x) in objByLayer[n - 1]])
            allLayerBuffObj = allLayerObj.buffer(self.u_width, join_style=2)
            resinRegions = [(aID, a) for (aID, a) in self._grids[n].iteritems()
                            if a.within(allLayerBuffObj)]
            for (i, (bID, b)) in enumerate(resinRegions):
                if b.objectType is None:
                    setattr(self._grids[n][bID], 'objectType', 'Resin')
                    self._grids[n][bID].layer = n
                    self.setAngle(self._grids[n][bID], passAngle)
                    resinRegions[i] = (bID, self._grids[n][bID])
                elif b.objectType == 'Undulation' and passAngle not in b.angle:
                    self.setAngle(self._grids[n][bID], passAngle, insert=True)

    def raiseTapes(self, passPath, passAngle):
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
        for i in range(len(passPath)):
            objID, obj = passPath[i]
            if obj.objectType is None:
                setattr(obj, 'objectType', 'Tape')
                self.setAngle(obj, passAngle)
            elif obj.objectType == 'Tape':
                # check layer above
                m = 2
                while m <= len(self._grids):
                    objGrid = self._grids[m]  # next grid layer
                    if objGrid[objID].objectType is None:
                        setattr(objGrid[objID], 'objectType', 'Tape')
                        objGrid[objID].layer = m
                        self.setAngle(objGrid[objID], passAngle)
                        passPath[i] = (objID, objGrid[objID])
                        break
                    else:
                        m += 1
        passPath = self.inbetween(passPath, passAngle, 'Tape')
        return passPath

    def inbetween(self, passPath, passAngle, objType):
        objByLayer = [
            [(objID, obj) for (objID, obj) in passPath
             if obj.layer == s and obj.objectType == objType]
            for s in range(1, len(self._grids) + 1)]
        for n in range(1, len(objByLayer)):
            # raise regions between objects in the same layer
            allObjInLayer = cascaded_union([x for (xID, x) in objByLayer[n]])
            try:
                allObjBuff = [obj.buffer(2 * self.u_width, join_style=2)
                              for obj in allObjInLayer.geoms]
                allObjConvexDiff = cascaded_union(
                    [a.intersection(b) for a, b in combinations(allObjBuff, 2)])
                for (zID, z) in passPath:
                    if z.within(allObjConvexDiff):
                        cell = self._grids[n + 1][zID]
                        cell.objectType = objType
                        cell.layer = n + 1
                        self.setAngle(cell, passAngle)
                        idx = [i for (i, el) in passPath].index(zID)
                        passPath[idx] = (zID, cell)
            except AttributeError:
                pass
        return passPath

    def createUndulations(self, passPath, passAngle):
        '''
        Create undulation regions
        '''
        objByLayer = [[(objID, obj) for (objID, obj)
                       in passPath if obj.layer == s] for s
                      in range(1, len(self._grids) + 1)]
        for n in range(1, len(objByLayer)):
            topLayer = cascaded_union([x for (xID, x) in objByLayer[n]])
            bottomLayer = objByLayer[n - 1]
            topLayerBuff = topLayer.buffer(2 * self.u_width + 1.0,
                                           join_style=2)
            bottomObjTouchingTopObj = [(zID, z) for (zID, z) in bottomLayer
                                       if z.within(topLayerBuff)]
            for (cellID, cell) in bottomObjTouchingTopObj:
                cellBottom = self._grids[n][cellID]
                cellTop = self._grids[n + 1][cellID]
                cellTop.objectType = 'Undulation'
                cellTop.layer = n + 1
                cellBottom.objectType = 'Undulation'
                if len(cellBottom.angle) != 2:
                    self.setAngle(cellBottom, passAngle)
                self.setAngle(cellTop, passAngle)
                passPath.append([cellID, cellTop])

        return passPath

    def defineConnectivity(self, passPath):
        connectivity = []
        for (objID, obj) in passPath:
            if objID in self.grids[int(obj.layer)].keys():
                ident = '{}-{}'.format(obj.layer, objID)
                setattr(obj, 'objectID', ident)
                connectivity.append(ident)
        return connectivity

    def setAngle(self, obj, angle, insert=False):
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
            if angle not in obj.angle:
                if insert:
                    obj.angle.insert(0, angle)
                else:
                    obj.angle.append(angle)


# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.


if __name__ == '__main__':
    from objectPlot import objPlot
    # import matplotlib.pyplot as plt

    tapeAngles = (0, 90)
    tapeWidths = (30.0, ) * len(tapeAngles)

    # initialize laminate
    interlacedLaminate = Interlaced(tapeAngles, tapeWidths)

    interlacedLaminate.makePass(passAngle=0)
    w = 30.0
    uw = 1.0
    crds = ([(-300.0, (w / 2) - uw - w), (300.0, (w / 2) - uw - w),
             (300.0, (-w / 2) + uw - w), (-300.0, (-w / 2) + uw - w)])
    interlacedLaminate.makePass(passAngle=0, passCoords=crds)
    interlacedLaminate.makePass(passAngle=90)
    crds = ([(-300.0, (w / 2) - uw - 2 * w), (300.0, (w / 2) - uw - 2 * w),
             (300.0, (-w / 2) + uw - 2 * w), (-300.0, (-w / 2) + uw - 2 * w)])
    interlacedLaminate.makePass(passAngle=0, passCoords=crds)
    grid = interlacedLaminate.returnTrimmedGrid()

    # for (idx, pol) in grid[1].iteritems():
    #     x, y = pol.exterior.xy
    #     plt.plot(x, y)
    # plt.show()

    objPlot(grid, tapeAngles, 'Tape')
    objPlot(grid, tapeAngles, 'Resin')
    objPlot(grid, tapeAngles, 'Undulation')


