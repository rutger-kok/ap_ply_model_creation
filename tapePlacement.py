import numpy as np 
import math 
from shapelyInterlacedGeometryCreation import Tape, Resin, Undulation 
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
'''
Define laminate using following 'gene'

(tapeAngles=(...), tapeWidths=(...), tapeSpacing=(...), t, r)
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

boundary = Polygon([(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])

# define tape width and thickness
w = 10.0
t = 0.2

# TODO code exception for 90 degree tapes
def laminateCreation(tapeAngles, tapeWidths, tapeSpacing, t=2.0, r=2.0):
    tapes = zip(tapeAngles, tapeWidths, tapeSpacing)
    for tape in tapes:
        a = tape[0]  # angle
        w = tape[1]  # width
        s = tape[2]  # spacing
        offset = (w*(1+s)+2*r)/math.cos(math.radians(a))
        maxOffset = (75.0+ (r+w/2.0)*math.cos(math.radians(a))
                + 50.0*math.tan(math.radians(a)))
        numberOffsets = int(maxOffset/offset)
        offsetList = [offset*k for k in range(-numberOffsets,numberOffsets+1)]

        for os in offsetList:
            tapeCoords = (
                    [(-1000.0, (-w/2.0)+os), (-1000.0, (w/2.0)+os),
                    (1000.0, (w/2.0)+os), (1000.0, (-w/2.0)+os)])
            resinCoords1 = (
                    [(-1000.0, (-w/2.0)+os-r), (-1000.0, (-w/2.0)+os),
                    (1000.0, (-w/2.0)+os), (1000.0, (-w/2.0)+os-r)])
            resinCoords2 = (
                    [(-1000.0, (w/2.0)+os), (-1000.0, (w/2.0)+os+r),
                    (1000.0, (w/2.0)+os+r), (1000.0, (w/2.0)+os)])
                                        
            # createdTape = Tape(coords=tapeCoords, angle=a)
            # createdResin1 = Resin(coords=resinCoords1, angle=a)
            # createdResin2 = Resin(coords=resinCoords2, angle=a)

laminate1 = laminateCreation(tapeAngles=(0,), tapeWidths=(10,), tapeSpacing=(1,))

# laminate1 = laminateCreation(tapeAngles=(0,90,45), tapeWidths=(10,10,10), tapeSpacing=(1,1,1))

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # Plotting of regions

    f1, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    f1.suptitle('Undulation per Layer')

    g = 1
    for i in range(3):
        for j in range(3):
            if list(Undulation.getinstances(g)):
                for k in range(len(list(Undulation.getinstances(g)))):
                    xi,yi = list(Undulation.getinstances(g))[k].exterior.xy
                    axes[i,j].plot(xi,yi)
                    axes[i,j].set_title('Layer: {}'.format(g))
                xb, yb = boundary.exterior.xy
                axes[i,j].plot(xb,yb)
            else:
                break
            g += 1

    f2, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    f2.suptitle('Tape per Layer')
    g = 1
    for i in range(3):
        for j in range(3):
            if list(Tape.getinstances(g)):
                for k in range(len(list(Tape.getinstances(g)))):
                    xi,yi = list(Tape.getinstances(g))[k].exterior.xy
                    axes[i,j].plot(xi,yi)
                    axes[i,j].set_title('Layer: {}'.format(g))
                xb, yb = boundary.exterior.xy
                axes[i,j].plot(xb,yb)
            else:
                break
            g += 1

    f3, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    f3.suptitle('Resin per Layer')

    g = 1
    for i in range(3):
        for j in range(3):
            if list(Resin.getinstances(g)):
                for k in range(len(list(Resin.getinstances(g)))):
                    xi,yi = list(Resin.getinstances(g))[k].exterior.xy
                    axes[i,j].plot(xi,yi)
                    axes[i,j].set_title('Layer: {}'.format(g))
                xb, yb = boundary.exterior.xy
                axes[i,j].plot(xb,yb)
            else:
                break
            g += 1
    plt.show(f1)
    plt.show(f2)
    plt.show(f3)