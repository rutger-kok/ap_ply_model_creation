from shapely.geometry import Polygon
from shapely import affinity
from shapely.geometry import Point
import numpy as np
from shapely.ops import cascaded_union
from itertools import count
from itertools import combinations
from shapely.geometry import shape, JOIN_STYLE

'''
This program is used to create interlaced laminate geometries with arbitrary
tape orientations. The coordinates of different regions are then used as 
inputs to Abaqus simulations. Three different regions are defined; tape 
regions, undulation regions, and pure resin regions.

Method:
The code defines three classes, Tape, Resin, and Undulation. The latter two are
subclasses of the Tape class, which is itself a subclass of the Polygon class
from the Shapely package. 

When a call is made to instantiate a Tape object, the program assigns the tape
by default to the first layer of the laminate. It then checks a list of 
existing Tape instances in that layer to determine whether the new Tape 
intersects with any ofthem. If an intersection is detected, then the region of
the new Tape that intersects is shifted up a layer (to layer 2). The program 
then checks whether the part of the tape that has been shifted up intersects 
with any Tape instances in that layer. This process repeats itself until the 
top of the laminate is reached (i.e. no intersections are detected).

The program will also check whether the Tape instance intersects with any 
existing Resin regions iin the current layer. If there are intersections, these
regions are identified as Undulation regions. The process of shifting the
intersecting regions up a layer is the same as for intersections with other
Tapes. 

The parts of the tapes that are not within the intersection region are placed
on the current layer. 

'''


# define outside boundary of the impact specimen
boundary = Polygon([(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])


# import matplotlib.pyplot as plt
# test = Polygon([(5.999982624777836, 13.071050436643302), (5.999993825921681, 13.071045796977769), (5.999993999999835, 6.000006222427927), (5.580271773394754, 5.999994632422353), (13.071046036769234, 5.9999923670399875), (12.07107307878631, 5.000005266920835), (5.99999665018369, 5.000006496179962), (5.999996970236527, 5.000006688012394), (5.999997364351884, 5.000006980308498), (5.999997727919422, 5.000007309827197), (5.999998057437783, 5.000007673395041), (5.999998349733522, 5.00000806751067), (5.999998601991668, 5.000008488378531), (5.999998811782834, 5.00000893194543), (5.999998977086614, 5.000009393939573), (5.999999096311042, 5.000009869911698), (5.999999168307918, 5.000010355277928), (5.999999192383874, 5.000010845363916), (5.999998301329779, 5.999987825139856), (5.999998277252977, 5.999988315225526), (5.999998205255316, 5.9999888005913595), (5.999998086030174, 5.99998927656302), (5.999997920725753, 5.999989738556639), (5.999997710934025, 5.999990182122967), (5.9999974586754, 5.999990602990217), (5.999997166379266, 5.999990997105209), (5.9999968368605945, 5.999991360672397), (5.999996473292829, 5.999991690190433), (5.999996079177328, 5.999991982485879), (5.999995658309636, 5.999992234743767), (5.999995214742942, 5.99999244453472), (5.999994752749034, 5.9999926098383325), (5.999994276777165, 5.999992729062643), (5.999993791411205, 5.999992801059455), (5.999993301325493, 5.999992825135401), (5.000017011540341, 5.999991968094749), (5.000016521454274, 5.999991944017925), (5.000016036088053, 5.999991872020164), (5.000015560116027, 5.9999917527948465), (5.000015098122073, 5.999991587490179), (5.000014654555454, 5.999991377698137), (5.000014233687964, 5.999991125439138), (5.000013839572795, 5.999990833142576), (5.000013476005497, 5.999990503623432), (5.000013146487428, 5.999990140055161), (5.000013038866356, 5.999989994944461), (5.00001144552027, 12.07107925738574), (5.999982624777836, 13.071050436643302)])
# test2 = Polygon([(-8.659311595910923, -4.999929775808399), (-2.886803213544549, 4.999910162687729), (5.0, 4.999923559819638), (5.0, 18.66025403784437), (5.999628046379085, 20.389711020306372), (6.0, 4.999925258496724), (8.659762545341952, 4.999148710014245), (4.123739840977385, -2.8574730782311635), (4.999773253753763, -1.340138698174518), (4.999660262690608, 4.999923557870636), (-2.886876749170912, 4.999782775346785), (-4.999838480402644, 1.339453901155823), (-4.99983847606975, -4.99976203561501), (2.8868964699636686, -4.999748637831682), (2.886803213544549, -4.999910162687729), (-5.0, -4.999923559819637), (-5.0, -18.66025403784437), (-5.998774257369033, -20.389853044684596), (-6.0, -20.3894022935068), (-6.0, -4.999925258496724), (-8.659311595910923, -4.999929775808399)])
# buffer1 = boundary.buffer(5, join_style=1)
# buffer2 = boundary.buffer(5, join_style=3)
# xb, yb = boundary.exterior.xy
# x1,y1 = buffer1.exterior.xy
# x2,y2 = buffer2.exterior.xy 
# plt.plot(xb,yb) 
# plt.plot(x1,y1, label='buffer 1')
# plt.plot(x2,y2, label='buffer 2')
# plt.legend()
# plt.show()

# define tape width and thickness
w = 10.0
t = 0.2

# define the centerpoint of the impact specimen
centerpoint = Point(0.0,0.0)

# defines the Tape class, which is a subclass of the Shapely Polygon class.
class Tape(Polygon):
    _ids = count(0)  # counts the number of instances created
    _instances = []  # initializes a list of Tape instances 
    def __init__(self, angle=0, coords=None , holes=None, layer=1,
            angleLabel=None, step1=True, check=True):
        self.createObject(angle, coords, holes, layer, angleLabel)
        if check:

            if step1:
                allObjInLayer = list(list(Tape.getinstances(self.layer)) + 
                        list(Undulation.getinstances(self.layer)) + 
                        list(Resin.getinstances(self.layer)))

                tapeIntAllCoordList, tapeSplitAll = self.checkIntersect(
                        allObjInLayer)
                if tapeIntAllCoordList:
                    for tapeIntAllCoords in tapeIntAllCoordList:
                        if self.clean(tapeIntAllCoords):
                            mergedTapeLayer = Tape(
                                coords=self.clean(tapeIntAllCoords),
                                layer=self.layer+1, angleLabel=self.angle)
                        else: continue
                    for tapeSplitAllCoords in self.splitObject(tapeSplitAll):
                        if self.clean(tapeSplitAllCoords):
                            tapeSplitAllObj = Tape(
                                coords=self.clean(tapeSplitAllCoords),
                                layer=self.layer, angleLabel=self.angle,
                                step1=False)     
                        else: continue           
                else:
                    if self.clean(self.coordinates):
                        Tape(coords=self.clean(self.coordinates), layer=self.layer, 
                            angleLabel = self.angle, step1=False)

            else:
                # check intersect with self (other Tape objects)
                tapeIntSelfCoordList, tapeSplitSelf = self.checkSelfIntersect()

                resinUndulationMerge = list(
                        list(Resin.getinstances(self.layer-1))+
                        list(Undulation.getinstances(self.layer-1)))
                tapeIntMergeCoordList, tapeSplitMerge = self.checkIntersect(
                        resinUndulationMerge)
                if tapeIntMergeCoordList:
                    for coordinateSet in tapeIntMergeCoordList:
                        if self.clean(coordinateSet):
                            newUndBottomLayer = Undulation(
                                coords=self.clean(coordinateSet),
                                layer=self.layer-1, angleLabel=self.angle)   
                            newUndTopLayer = Undulation(
                                coords=self.clean(coordinateSet),
                                layer=self.layer, angleLabel=self.angle)
                        else: continue

                        for resinObject in Resin.getinstances(self.layer-1):
                            diffObject = Polygon(coordinateSet)
                            if diffObject.intersection(resinObject):
                                Resin._instances.remove(resinObject)
                                resinSplitObjectList = resinObject.difference(
                                        diffObject)#.buffer(0.01))
                                resinSplitCoords = self.splitObject(
                                        resinSplitObjectList)       
                                # create split resin regions
                                for coordinateSet2 in resinSplitCoords:
                                    if self.clean(coordinateSet2):
                                        splitResin = Resin(
                                            coords=self.clean(coordinateSet2),
                                            layer=self.layer-1,
                                            angleLabel=self.angle, check=False)
                                    else: continue

                # define region of Tape not intersecting with any other object 
                # (not with other Tape, Undulation, or Resin objects)
                totalSplitObjList = None
                for item in (tapeSplitMerge, tapeSplitSelf):
                    if item and not item.is_empty:
                        if totalSplitObjList:
                            totalSplitObjList = totalSplitObjList.intersection(
                                    item)#.buffer(-1*10**-4))
                        else:
                            totalSplitObjList = item

                # create split tape object (for non-intersecting Tape regions)
                if totalSplitObjList:
                    for tapeSplitCoords in self.splitObject(totalSplitObjList):
                        if self.clean(tapeSplitCoords):
                            tapeSplitObj = self.newInstance(
                                tempCoords=self.clean(tapeSplitCoords),
                                tempLayer=self.layer,
                                tempAngleLabel=self.angle)
                        else: continue

                # if the current Tape doesnt intersect with any other object it
                # is added to the instance list
                if not tapeIntSelfCoordList and not tapeIntMergeCoordList:
                    for pgon in list(self.getinstances(self.layer)):
                        if self.almost_equals(pgon, decimal=3):
                            break
                    else:
                        if self.area > 0.001:
                            self._instances.append(self)
        else:
            for pgon in list(self.getinstances(self.layer)):
                if self.almost_equals(pgon, decimal=3):
                    break
            else:
                if self.area > 0.001:
                    self._instances.append(self)

    def clean(self, coordinates):
        poly = Polygon(coordinates)
        tol = 0.001 # distance
        cf = 1.3  # cofactor
        cleaned = poly.buffer(-tol).buffer(tol*cf).intersection(poly).simplify(tol)
        try:
            new_coords = list(zip(cleaned.exterior.xy[0],cleaned.exterior.xy[1]))
            return new_coords
        except AttributeError:
            print coordinates
            print cleaned
            print '__________________'
            return None
        

    def checkSelfIntersect(self):
        # first check intersection with own type e.g. for tapes check 
        # intersect with other tapes
        intersectSelfList, splitSelfList = self.checkIntersect(
                list(self.getinstances(self.layer)))
        if intersectSelfList:
            for intersectSelfCoords in intersectSelfList:
                if self.clean(intersectSelfCoords):
                    intersectSelfObj = self.newInstance(
                        tempCoords=self.clean(intersectSelfCoords),
                        tempLayer=self.layer+1, tempAngleLabel=self.angle)
                else: continue
        return intersectSelfList, splitSelfList

    def createObject(self, angle=0, coords=None , holes=None, layer=1,
            angleLabel=None):
        if coords == None:
            coords = (
                    [(-100.0, -w/2), (100.0, -w/2), (100.0, w/2), (-100.0, w/2)]
                    )
        a = Polygon(coords, holes)
        b = affinity.rotate(a, angle, origin=centerpoint)
        c = b.intersection(boundary) # trim tape using boundary
        x0r = map(lambda x: round(x, 5), c.exterior.xy[0]) 
        y0r = map(lambda x: round(x, 5), c.exterior.xy[1])
        new_coords = zip(x0r,y0r)
        Polygon.__init__(self,new_coords, holes) # initialize rotated tape

        # object parameters
        self.coordinates = new_coords  
        self.layer = layer
        self.angle = angle
        self.objNumber = next(self._ids)
        
        # need to define another variable for angle that does not rotate the 
        # object this is necessary to define material orientations in abaqus 
        if angleLabel == None:  
            self.angleLabel = self.angle
        else:
            self.angleLabel = angleLabel

    # this class method is a generator function used to return the instances of
    # class in a specific requested layer
    @classmethod
    def getinstances(cls, layer):
        dead = []
        for obj in cls._instances:
            if obj is not None:
                if obj.layer == layer:
                    yield obj
            else:
                dead.append(obj)
        cls._instances = [x for x in cls._instances if x not in dead]

    # this class method is used to create new instances of whatever class it is
    # called from
    @classmethod  
    def newInstance(cls, tempCoords, tempLayer, tempAngleLabel):
        return cls(
                coords=tempCoords, layer=tempLayer, angleLabel=tempAngleLabel)

    # this method checks whether the current object intersects with any of the
    # objects in "objectList"
    def checkIntersect(self, objectList, buffer=True):
        differenceObj = None
        intersectCoords = []
        if objectList:
            if buffer:  # applies negative buffer if requested
                # note the objects in the objectList are merged before the 
                # intersection check
                mergedObj = cascaded_union(objectList).buffer(-1*10**-5, join_style=2)
            else:
                mergedObj = cascaded_union(objectList)
            intersectObj = self.intersection(mergedObj)
            if intersectObj:
                # different behaviour depending on whether the intersection
                # creates multiple intersecting regions or only one
                if intersectObj.geom_type == 'MultiPolygon':
                    for plygn in intersectObj:
                        intersectCoords.append(
                                zip(plygn.exterior.xy[0], plygn.exterior.xy[1]))
                        intersectObjBuffered = intersectObj.buffer(1*10**-5, join_style=2)
                        differenceObj = self.difference(intersectObjBuffered)
                elif intersectObj.geom_type == 'Polygon':
                    intersectCoords.append(
                            zip(intersectObj.exterior.xy[0], 
                            intersectObj.exterior.xy[1]))
                    intersectObjBuffered = intersectObj.buffer(1*10**-5, join_style=2)
                    differenceObj = self.difference(intersectObjBuffered)
                else:
                    print 'Not a polygon or multipolygon'     
        return intersectCoords, differenceObj

    # this method is used to define the coordinates of the objects
    # "left behind" when the current object intersects with another object
    def splitObject(self, differenceObject):
        splitCoords = []
        if differenceObject:
            if differenceObject.geom_type == 'MultiPolygon':
                for plg in differenceObject:
                    splitCoords.append(
                            zip(plg.exterior.xy[0], plg.exterior.xy[1]))
            elif differenceObject.geom_type == 'Polygon':
                splitCoords.append(
                        zip(differenceObject.exterior.xy[0],
                        differenceObject.exterior.xy[1]))
            else:
                print 'Not a polygon or multipolygon'
        return splitCoords

    def __repr__(self):
        return 'Tape: Layer: {}, Number: {}'.format(self.layer, self.objNumber)
    def __str__(self):
        return 'Tape: Layer: {}, Number: {}'.format(self.layer, self.objNumber)               

# the Resin class is a subclass of the Tape class. Aside from object creation, 
# the class defines its own __init__ method as its intersection behaviour is
# different (resin intersecting tape -> no undulation). 
class Resin(Tape):
    _ids = count(0)
    _instances = []

    def __init__(self, angle=0, coords=None , holes=None, layer=1,
            angleLabel=None, check=True):
        self.createObject(angle, coords, holes, layer, angleLabel)

        if check:
            resinIntSelfCoordList, resinSplitSelf = self.checkSelfIntersect()

            # check intersect with other e.g. for tapes check intersect with resin 
            resinIntTapeCoordList, resinSplitTape = self.checkIntersect(
                    list(Tape.getinstances(self.layer)))

            if resinIntTapeCoordList:
                for resinIntTapeCoords in resinIntTapeCoordList:
                    if self.clean(resinIntTapeCoords):
                        resinIntTapeObj = self.newInstance(
                            tempCoords=self.clean(resinIntTapeCoords),
                            tempLayer=self.layer+1, tempAngleLabel=self.angle)
                    else: continue

            resinIntUndCoordList, resinSplitUnd = self.checkIntersect(
                    list(Undulation.getinstances(self.layer)))                            
            if resinIntUndCoordList:
                for resinIntUndCoords in resinIntUndCoordList:
                    if self.clean(resinIntUndCoords):
                        resinIntUndObj = self.newInstance(
                            tempCoords=self.clean(resinIntUndCoords),
                            tempLayer=self.layer+1, tempAngleLabel=self.angle)
                    else: continue

            totalSplitObjList = None
            for item in (resinSplitSelf, resinSplitTape, resinSplitUnd):
                if item and not item.is_empty:
                    if totalSplitObjList:
                        totalSplitObjList = totalSplitObjList.intersection(
                                item)#.buffer(-1*10**-3, join_style=2))
                    else:
                        totalSplitObjList = item

            if totalSplitObjList:
                for resinSplitCoords in self.splitObject(totalSplitObjList):
                    if self.clean(resinSplitCoords):
                        resinSplitObj = self.newInstance(
                            tempCoords=self.clean(resinSplitCoords),
                            tempLayer=self.layer, tempAngleLabel=self.angle)
                    else: continue                      

            if not resinIntSelfCoordList and not resinIntUndCoordList and not resinIntTapeCoordList:
                for pgon in list(self.getinstances(self.layer)):
                    if self.almost_equals(pgon, decimal=3):
                        break
                else:
                    if self.area > 0.001:
                        self._instances.append(self)
        else:
            if self.area > 0.001:
                self._instances.append(self)

    def __repr__(self):
        return 'Resin: Layer: {}, Number: {}'.format(self.layer, self.objNumber)
    def __str__(self):
        return 'Resin: Layer: {}, Number: {}'.format(self.layer, self.objNumber)   

# the Undulation class is a subclass of the Tape class. The Undulation class
# defines its own __init__ method which defines how the undulation regions
# are recorded.
class Undulation(Tape):
    _ids = count(0)
    _instances = []
    def __init__(self, angle=0, coords=None , holes=None, layer=1,
            angleLabel=None, check=True):
        self.createObject(angle, coords, holes, layer, angleLabel)
        if check:
            for prevUnd in list(self.getinstances(self.layer)):
                if self.intersection(prevUnd):
                    prevSplitUndObj = None
                    newSplitUndulationObj = None
                    intersectCoords = []
                    selfBuff = (self.buffer(1*10**-10,join_style=2))
                    prevUndBuff = prevUnd.buffer(-1*10**-10, join_style=2)
                    intersectObj = selfBuff.intersection(prevUndBuff)
                    if intersectObj:
                        if intersectObj.geom_type == 'MultiPolygon':
                            for plygn in intersectObj:
                                intersectCoords.append(
                                        zip(plygn.exterior.xy[0],
                                        plygn.exterior.xy[1]))
                                # intersectObjBuffered = intersectObj.buffer(
                                #         1*10**-5, join_style=2)
                                prevSplitUndObj = (prevUndBuff.difference(
                                        intersectObj))
                                newSplitUndulationObj = selfBuff.difference(
                                        intersectObj)
                        elif intersectObj.geom_type == 'Polygon':
                            intersectCoords.append(
                                    zip(intersectObj.exterior.xy[0],
                                    intersectObj.exterior.xy[1]))
                            # intersectObjBuffered = intersectObj.buffer(
                            #         1*10**-5, join_style=2)
                            prevSplitUndObj = (prevUndBuff.difference(
                                    intersectObj))
                            newSplitUndulationObj = selfBuff.difference(
                                    intersectObj)

                        prevSplitUndulationCoords = self.splitObject(
                                prevSplitUndObj)
                        newSplitUndulationCoords = self.splitObject(
                                newSplitUndulationObj)
                        for coordinateSet in intersectCoords:
                            if self.clean(coordinateSet):
                                newUndulation = Undulation(
                                    coords=self.clean(coordinateSet),
                                    layer=self.layer,
                                    angleLabel=(self.angle,
                                        prevUnd.angle), check=False)  
                            else: continue           
                        for coordinateSet2 in prevSplitUndulationCoords:
                            if self.clean(coordinateSet2):
                                prevSplitUndulation = Undulation(
                                    coords=self.clean(coordinateSet2),
                                    layer=self.layer,
                                    angleLabel=prevUnd.angle,
                                    check=False)             
                            else: continue
                        for coordinateSet3 in newSplitUndulationCoords:
                            if self.clean(coordinateSet3):
                                newSplitUndulation = Undulation(
                                    coords=self.clean(coordinateSet3), 
                                    layer=self.layer, 
                                    angleLabel=self.angle,
                                    check=False)            
                            else: continue
                        
                        self._instances.remove(prevUnd)
                else:
                    continue
            else:
                if self.clean(self.coordinates):
                    newUndulation = Undulation(
                        coords=self.clean(self.coordinates), layer=self.layer,
                        angleLabel=self.angle, check=False)

        else:
            for pgon in list(self.getinstances(self.layer)):
                if self.almost_equals(pgon, decimal=1):
                    continue
            else:
                if self.area > 0.001:
                    self._instances.append(self)

    def __repr__(self):
        return 'Undulation: Layer: {}, Number: {}'.format(
                self.layer, self.objNumber)
    def __str__(self):
        return 'Undulation: Layer: {}, Number: {}'.format(
                self.layer, self.objNumber)   

tape1 = Tape(angle=0)
testResin2 = Resin(
        coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)])
testResin2 = Resin(
        coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)])

tape3 = Tape(angle=90)
testResin3 = Resin(
        coords=[(-6.0, 75.0), (-5.0, 75.0), (-5.0, -75.0), (-6.0, -75.0)])
testResin4 = Resin(
        coords=[(5.0, 75.0), (6.0, 75.0), (6.0, -75.0), (5.0, -75.0)])

tape4 = Tape(coords=[(6.0, 75.0), (16.0, 75.0), (16.0, -75.0), (6.0, -75.0)])
# tape2 = Tape(angle=60)
# tape3 = Tape(angle=-60)

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
    
    # print Tape._instances
    # print Resin._instances
    # print Undulation._instances

    # f4, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    # f4.suptitle('Undulation per Layer')
    # g = 0
    # for i in range(3):
    #     for j in range(3):
    #         try:
    #             xi,yi = list(Tape.getinstances(2))[g].exterior.xy
    #             axes[i,j].plot(xi,yi)
    #             axes[i,j].set_title('Layer: {}'.format(g))
    #             xb, yb = boundary.exterior.xy
    #             axes[i,j].plot(xb,yb)
    #             g += 1
    #         except IndexError:
    #             break

    # plt.show(f4)


