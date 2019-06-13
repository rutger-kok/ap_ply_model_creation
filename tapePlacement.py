import numpy as np 
import math 
import shapelyInterlacedGeometryCreation as sigc
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

# define specimen boundaries
boundary = Polygon([(-50.1, -75.1), (50.1, -75.1), (50.1, 75.1), (-50.1, 75.1)])

# define tape width and thickness
w = 10.0
t = 0.2

def laminateCreation(tapeAngles, tapeWidths, tapeSpacing, t=2.0, r=1.0):
    tapes = zip(tapeAngles, tapeWidths, tapeSpacing)

    for tape in tapes:
        a = tape[0]  # angle
        w = tape[1]  # width
        s = tape[2]  # spacing
        for numShift in range(s+1):
            if a == 90:
                offset = (w*(1+s)+2*r)
                maxOffset = 50.0+r+w/2.0
                numberOffsets = int(maxOffset/offset)
                offsetList = [offset*k for k in range(-numberOffsets-numShift,numberOffsets+1-numShift)]
                shift = numShift*offset
                for os in offsetList:
                    tapeCoords = (
                            [((-w/2.0)+os+shift, 100.0),
                            ((w/2.0)+os+shift, 100.0),
                            ((w/2.0)+os+shift, -100.0),
                            ((-w/2.0)+os+shift, -100.0)])
                    resinCoords1 = (
                            [((-w/2.0)+os-r+shift, 100.0),
                            ((-w/2.0)+os+shift, 100.0),
                            ((-w/2.0)+os+shift, -100.0),
                            ((-w/2.0)+os-r+shift, -100.0)])
                    resinCoords2 = (
                            [((w/2.0)+os+shift, 100.0),
                            ((w/2.0)+os+r+shift, 100.0),
                            ((w/2.0)+os+r+shift, -100.0),
                            ((w/2.0)+os+shift, -100.0)])
                                            
                    createdTape = sigc.machinePass(coords=tapeCoords)
                    createdResin1 = sigc.Resin(coords=resinCoords1)
                    createdResin2 = sigc.Resin(coords=resinCoords2)
            else:
                offset = (w*(1+s)+2*r)/math.cos(math.radians(a))
                maxOffset = (75.0+ (r+w/2.0)*math.cos(math.radians(a))
                        + 50.0*math.tan(math.radians(a)))
                numberOffsets = int(maxOffset/offset)
                offsetList = [offset*k for k in range(-numberOffsets-numShift,numberOffsets+1-numShift)]
                shift = numShift*offset
                for os in offsetList:
                    tapeCoords = (
                            [(-1000.0, (w/2.0)+os+shift),
                            (1000.0, (w/2.0)+os+shift),
                            (1000.0, (-w/2.0)+os+shift),
                            (-1000.0, (-w/2.0)+os+shift)])
                    resinCoords1 = (
                            [(-1000.0, (-w/2.0)+os-r+shift),
                            (-1000.0, (-w/2.0)+os+shift),
                            (1000.0, (-w/2.0)+os+shift),
                            (1000.0, (-w/2.0)+os-r+shift)])
                    resinCoords2 = (
                            [(-1000.0, (w/2.0)+os+shift),
                            (-1000.0, (w/2.0)+os+r+shift),
                            (1000.0, (w/2.0)+os+r+shift),
                            (1000.0, (w/2.0)+os+shift)])
                                            
                    createdTape = sigc.machinePass(coords=tapeCoords, angle=a)
                    rotPoint = createdTape.tapeCtd
                    createdResin1 = sigc.Resin(coords=resinCoords1, angle=a, rPoint=rotPoint)
                    createdResin2 = sigc.Resin(coords=resinCoords2, angle=a, rPoint=rotPoint)

laminate1 = laminateCreation(tapeAngles=(0,90), tapeWidths=(10,10), tapeSpacing=(1,1))

# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.

if __name__ == '__main__':
    import matplotlib.pyplot as plt
# # Plotting of regions
# from shapely.geometry import Polygon
# boundary = Polygon([(-50.1, -75.1), (50.1, -75.1), (50.1, 75.1), (-50.1, 75.1)])

# # x = [941.4027215025367, -937.9825200692801, -938.3245402126057, 941.0607013592111, 941.4027215025367]
# # y = [241.54568069890706, -442.49460595243033, -441.5549133316444, 242.485373319693, 241.54568069890706]
# x1,y1 = boundary.exterior.xy

# # plt.plot(x,y)
# plt.plot(x1,y1)
# # plt.show()
# b = Polygon([(-39.00000005, 29.00000003), (-39.00000005, 29.00000002), (-39.00000005, 38.99999998), (-39.00000005, 38.99999997), (-39.00000005, 38.99999998), (-39.00000007, 38.99999998), (-39.00000007, 39.00000005), (-39.00000001, 39.00000005), (-39.00000001, 39.00000004), (-39.0, 39.00000004), (-39.0, 39.00000004), (-29.0, 39.00000004), (-29.0, 39.00000004), (-29.0, 39.00000004), (-28.99999999, 39.00000004), (-28.99999999, 39.00000005), (-28.99999995, 39.00000005), (-28.99999995, 39.00000004), (-28.99999993, 39.00000004), (-28.99999993, 38.99999998), (-28.99999995, 38.99999998), (-28.99999995, 38.99999997), (-28.99999995, 38.99999998), (-28.99999995, 29.00000002), (-28.99999995, 29.00000003), (-28.99999995, 29.00000002), (-28.99999993, 29.00000002), (-28.99999993, 28.99999996), (-28.99999995, 28.99999996), (-28.99999995, 28.99999995), (-28.99999999, 28.99999995), (-28.99999999, 28.99999996), (-29.0, 28.99999996), (-29.0, 28.99999996), (-29.0, 28.99999996), (-39.0, 28.99999996), (-39.0, 28.99999996), (-39.0, 28.99999996), (-39.00000001, 28.99999996), (-39.00000001, 28.99999995), (-39.00000005, 28.99999995), (-39.00000005, 28.99999996), (-39.00000007, 28.99999996), (-39.00000007, 29.00000002), (-39.00000005, 29.00000002), (-39.00000005, 29.00000003)])
# xb,yb = b.exterior.xy 
# plt.plot(xb,yb)
# plt.show()
# print Polygon(zip(x,y)).centroid

    f1, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    f1.suptitle('Undulation per Layer')

    g = 1
    for i in range(3):
        for j in range(3):
            if list(sigc.Undulation.getinstances(g)):
                for k in range(len(list(sigc.Undulation.getinstances(g)))):
                    xi,yi = list(sigc.Undulation.getinstances(g))[k].exterior.xy
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
            if list(sigc.Tape.getinstances(g)):
                for k in range(len(list(sigc.Tape.getinstances(g)))):
                    xi,yi = list(sigc.Tape.getinstances(g))[k].exterior.xy
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
            if list(sigc.Resin.getinstances(g)):
                for k in range(len(list(sigc.Resin.getinstances(g)))):
                    xi,yi = list(sigc.Resin.getinstances(g))[k].exterior.xy
                    axes[i,j].plot(xi,yi)
                    axes[i,j].set_title('Layer: {}'.format(g))
                xb, yb = boundary.exterior.xy
                axes[i,j].plot(xb,yb)
            else:
                break
            g += 1

    # print sigc.Resin._instances
    
    plt.show(f1)
    plt.show(f2)
    plt.show(f3)

