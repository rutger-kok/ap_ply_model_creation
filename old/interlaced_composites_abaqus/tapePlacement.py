import math
import sigc

'''
Define laminate using following parameters

(tapeAngles=(...), tapeWidths=(...), tapeSpacing, r)
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


def laminateCreation(grid, tapeAngles, tapeWidths, tapeSpacing,
                     undulationWidth=1.0):
    paths = []
    tapes = zip(tapeAngles, tapeWidths)
    s = tapeSpacing  # spacing
    r = undulationWidth
    for numShift in range(s+1):
        for tape in tapes:
            a = tape[0]  # angle
            w = tape[1]  # width
            if a == 90:
                offset = w
                spacedOffset = (1+s)*offset
                maxOffset = 50.0+w/2.0
                numberOffsets = int(maxOffset/spacedOffset)
                offsetList = [spacedOffset*k for k
                              in range(-numberOffsets-numShift,
                                       numberOffsets+3-numShift)]
                shift = numShift*offset
                for os in offsetList:
                    tapeCoords = (
                                  [(-100.0+os+shift, (w/2.0)-r),
                                   (100.0+os+shift, (w/2.0)-r),
                                   (100.0+os+shift, (-w/2.0)+r),
                                   (-100.0+os+shift, (-w/2.0)+r)])
                    createdTape = sigc.machinePass(
                            grid, coords=tapeCoords, angle=90,
                            undulationWidth=r).connectivity
                    paths.append(createdTape)
            else:
                offset = w/math.cos(math.radians(a))
                spacedOffset = ((1+s)*w)/math.cos(math.radians(a))
                maxOffset = ((75.0+w/2.0)*math.cos(math.radians(a))
                             + 50.0*math.tan(math.radians(a)))
                numberOffsets = int(maxOffset/spacedOffset)
                offsetList = [spacedOffset*k for k
                              in range(-numberOffsets-numShift-5,
                                       numberOffsets-numShift+5)]
                shift = numShift*offset
                for os in offsetList:
                    tapeCoords = (
                            [(-150.0, ((w/2.0)-r)+os+shift),
                             (150.0, ((w/2.0)-r)+os+shift),
                             (150.0, ((-w/2.0)+r)+os+shift),
                             (-150.0, ((-w/2.0)+r)+os+shift)])
                    createdTape = sigc.machinePass(
                            grid, coords=tapeCoords, angle=a,
                            undulationWidth=r).connectivity
                    paths.append(createdTape)
    return paths

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.


if __name__ == '__main__':
    from objectPlot import objPlot

    # define tape width and thickness
    tapeAng = (0, 90)
    tapeW = (25,)*len(tapeAng)
    cpt = 0.2
    undulationRatio = 0.2
    rw = cpt/undulationRatio

    grid1 = sigc.createGrids(tapeAngles=tapeAng, tapeWidths=tapeW,
                             undulationWidth=rw)

    tapePaths = laminateCreation(grid1, tapeAngles=tapeAng,
                                 tapeWidths=tapeW, tapeSpacing=1,
                                 undulationWidth=rw)

    objPlot(grid1, tapeAng, 'Tape')
    objPlot(grid1, tapeAng, 'Resin')
    objPlot(grid1, tapeAng, 'Undulation')
