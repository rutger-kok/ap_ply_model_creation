import numpy as np 
import math 
import sigc3 as sigc
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

'''
Define laminate using following 'gene'

(tapeAngles=(...), tapeWidths=(...), tapeSpacing, r)
where:
tapeAngles = tuple of tape angles (floats) corresponding to the number of tapes
note that the tapes will be 'placed' in the sequence that they are defined!
tapeWidths = tuple of tape widths corresponding to each tape defined in the 
tapeAngles parameter (i.e. first tape width corresponds to the tape with the 
first tape angle). Note that 0 corresponds to NO gap, 1= 1 gap between tapes.
tapeSpacing = tuple of integers (1 per tape) defining the number of tape gaps 
to leave between tapes (as above, each spacing term corresponds to the tape 
with the same index)
t = tape thickness, this must be the same for each tape 
r = resin width assumed to be the same for each tape
'''

def laminateCreation(grid, tapeAngles, tapeWidths, tapeSpacing, r=1.0):
    tapes = zip(tapeAngles, tapeWidths)
    s = tapeSpacing  # spacing
    for numShift in range(s+1):
        for tape in tapes:
            a = tape[0]  # angle
            w = tape[1]  # width
            if a == 90:
                offset = w+2*r
                spacedOffset = (1+s)*offset
                maxOffset = 50.0+r+w/2.0
                numberOffsets = int(maxOffset/spacedOffset)
                offsetList = [spacedOffset*k for k 
                    in range(-numberOffsets-numShift,numberOffsets+2-numShift)]
                shift = numShift*offset
                for os in offsetList:
                    tapeCoords = (
                            [(-100.0+os+shift, (w/2.0)),
                            (100.0+os+shift, (w/2.0)),
                            (100.0+os+shift, (-w/2.0)),
                            (-100.0+os+shift, (-w/2.0))])   
                    createdTape = sigc.machinePass(
                            grid, coords=tapeCoords, angle=90)
            else:
                offset = ((w+2*r))/math.cos(math.radians(a))
                spacedOffset = ((1+s)*(w+2*r))/math.cos(math.radians(a))
                maxOffset = (75.0+ (r+w/2.0)*math.cos(math.radians(a))
                        + 50.0*math.tan(math.radians(a)))
                numberOffsets = int(maxOffset/spacedOffset)
                offsetList = [spacedOffset*k for k 
                    in range(-numberOffsets-numShift,numberOffsets+1-numShift)]
                shift = numShift*offset
                for os in offsetList:
                    tapeCoords = (
                            [(-150.0, (w/2.0)+os+shift),
                            (150.0, (w/2.0)+os+shift),
                            (150.0, (-w/2.0)+os+shift),
                            (-150.0, (-w/2.0)+os+shift)])    
                    createdTape = sigc.machinePass(
                            grid, coords=tapeCoords, angle=a)


# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.

if __name__ == '__main__':
    # define tape width and thickness
    grid1 = sigc.createGrids(tapeAngles=(0,90), tapeWidths=(10,10))
    laminateCreation(grid1, tapeAngles=(0,90), tapeWidths=(10,10), tapeSpacing=2)
    sigc.objPlot(grid1,'Tape')

