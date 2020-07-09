from math import cos, tan, radians
import sigc as sigc
from shapely.geometry import Polygon

'''
Define laminate using following parameters

(tapeAngles=(...), tapeWidths=(...), tapeSpacing, uw)
where:
tapeAngles = tuple of tape angles (floats) corresponding to the number of tapes
note that the tapes will be 'placed' in the sequence that they are defined!
tapeWidths = tuple of tape widths corresponding to each tape defined in the
tapeAngles parameter (i.e. first tape width corresponds to the tape with the
first tape angle).
tapeSpacing = tuple of integers (1 per tape) defining the number of tape gaps
to leave between tapes (as above, each spacing term corresponds to the tape
with the same index) Note that 0 corresponds to NO gap, 1->1 gap between tapes.

'''


def laminateCreation(tapeAngles, tapeWidths, tapeSpacing, specimenSize,
                     xMax=50.0, yMax=75.0, undulationWidth=1.0):
    # create interlaced laminate instance
    laminate = sigc.Interlaced(tapeAngles, tapeWidths, undulationWidth,
                               specimen=specimenSize)
    paths = []
    tapeProps = zip(tapeAngles, tapeWidths)
    s = tapeSpacing  # spacing
    uw = undulationWidth
    for numShift in range(s + 1):
        for (a, w) in tapeProps:
            if a == 90:
                offset = w
                spacedOffset = (1 + s) * offset
                maxOffset = xMax + w / 2.0
                numOffsets = int(maxOffset / spacedOffset)
                offsetList = [spacedOffset * k for k
                              in range(-numOffsets - numShift - 10,
                                       numOffsets - numShift + 10)]
                shift = numShift * offset
                for os in offsetList:
                    tapeCoords = ([(-300.0 + os + shift, (w / 2.0) - uw),
                                   (300.0 + os + shift, (w / 2.0) - uw),
                                   (300.0 + os + shift, (-w / 2.0) + uw),
                                   (-300.0 + os + shift, (-w / 2.0) + uw)])
                    createdTape = laminate.makePass(passCoords=tapeCoords,
                                                    passAngle=90)
                    paths.append(createdTape)
            else:
                offset = w / cos(radians(a))
                spacedOffset = ((1 + s) * w) / cos(radians(a))
                maxOffset = ((yMax + w / 2.0) * cos(radians(a))
                              +  xMax * tan(radians(a)))
                numOffsets = int(maxOffset / spacedOffset)
                offsetList = [spacedOffset * k for k
                              in range(-numOffsets - numShift - 10,
                                       numOffsets - numShift + 10)]
                shift = numShift * offset
                for os in offsetList:
                    tapeCoords = ([(-300.0, ((w / 2.0) - uw) + os + shift),
                                   (300.0, ((w / 2.0) - uw) + os + shift),
                                   (300.0, ((-w / 2.0) + uw) + os + shift),
                                   (-300.0, ((-w / 2.0) + uw) + os + shift)])
                    createdTape = laminate.makePass(passCoords=tapeCoords,
                                                    passAngle=a)
                    paths.append(createdTape)
    trimmedGrid = laminate.returnTrimmedGrid()    
    return paths, trimmedGrid

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.


if __name__ == '__main__':
    from objectPlot import objPlot
    # import matplotlib.pyplot as plt

    # define tape width and thickness
    tapeAngles = (0, 90)
    tapeWidths = (20,) * len(tapeAngles)
    tapeSpacing = 1
    cpt = 0.18
    undulationRatio = 0.18
    rw = cpt / undulationRatio

    xMin = yMin = -(tapeWidths[0] / 2.0)
    xMax = yMax = xMin + (tapeSpacing + 1) * (tapeWidths[0])
    numLayers = len(tapeAngles) * 2.0  # symmetric
    laminateThickness = numLayers * cpt
    zMax = laminateThickness / 2.0
    zMin = -zMax
    dimensions = [xMin, yMin, zMin, xMax, yMax, zMax]

    rvePoly = Polygon([(xMin, yMin), (xMax, yMin), (xMax, yMax), (xMin, yMax)])

    tapePaths, grid = laminateCreation(
        tapeAngles, tapeWidths, tapeSpacing, rvePoly,
        undulationWidth=rw, xMax=xMax, yMax=yMax)

    # for (key,poly) in grid[1].iteritems():
    #     x,y = poly.exterior.xy
    #     plt.plot(x,y)
    # plt.grid(True)
    # plt.show()

    objPlot(grid, tapeAngles, 'Tape', sample=rvePoly)
    objPlot(grid, tapeAngles, 'Resin', sample=rvePoly)
    objPlot(grid, tapeAngles, 'Undulation', sample=rvePoly)
