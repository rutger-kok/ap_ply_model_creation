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
# import matplotlib.pyplot as plt

# test = Polygon([(5.0, 5.0008), (5.0, 5.0008), (5.0, 5.0), (5.0, 5.0), (5.0, -1.33974), (2.88675, -5.0), (-5.0, -5.0), (-5.0, -5.0008), (-5.0, -5.0008), (-5.0, -5.0), (-5.0, -5.0), (-5.0, 1.33974), (-2.88675, 5.0), (5.0, 5.0), (5.0, 5.0008)])
# # define outside boundary of the impact specimen
boundary = Polygon([(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])
# xb, yb = boundary.exterior.xy
# x1,y1 = test.exterior.xy 
# plt.plot(xb,yb)
# plt.plot(x1,y1)
# plt.show()

# # define tape width and thickness
w = 10.0
t = 0.2

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
                            diffObject = Polygon(coordinateSet).buffer(0, join_style=2)
                            if diffObject.intersection(resinObject):
                                Resin._instances.remove(resinObject)
                                resinSplitObjectList = resinObject.difference(
                                        diffObject)
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
                                    item)#.buffer(-1*10**-6))
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
        tol = 1*10**-4 # distance
        cf = 1.5  # cofactor
        cleaned = poly.buffer(-tol,join_style=2).buffer(
                tol*cf, join_style=2).intersection(poly)#.simplify(1*10**-10)
        try:
            new_coords = list(zip(cleaned.exterior.xy[0],cleaned.exterior.xy[1]))
            return new_coords
        except AttributeError:
            # print coordinates
            # print cleaned
            # print '__________________'
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
        b = affinity.rotate(a, angle, origin=a.centroid)
        c = b.intersection(boundary) # trim tape using boundary
        # x0r = map(lambda x: round(x, 6), c.exterior.xy[0]) 
        # y0r = map(lambda x: round(x, 6), c.exterior.xy[1])
        # new_coords = zip(x0r,y0r)
        new_coords = zip(c.exterior.xy[0],c.exterior.xy[1])
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
    def checkIntersect(self, objectList):
        differenceObj = None
        intersectCoords = []
        if objectList:
            mergedObj = cascaded_union(objectList).buffer(-1*10**-6, join_style=2)
            selfBuffer = self.buffer(0, join_style=2)
            intersectObj = selfBuffer.intersection(mergedObj)
            if intersectObj:
                # different behaviour depending on whether the intersection
                # creates multiple intersecting regions or only one
                if intersectObj.geom_type == 'MultiPolygon':
                    for plygn in intersectObj:
                        intersectCoords.append(
                                zip(plygn.exterior.xy[0], plygn.exterior.xy[1]))
                        differenceObj = self.difference(intersectObj)
                elif intersectObj.geom_type == 'Polygon':
                    intersectCoords.append(
                            zip(intersectObj.exterior.xy[0], 
                            intersectObj.exterior.xy[1]))
                    differenceObj = self.difference(intersectObj)
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
        return 'Tape: Layer: {}, Number: {}, Coordinates: {}'.format(self.layer, self.objNumber, self.coordinates)
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
                    selfBuff = (self.buffer(0,join_style=2))
                    prevUndBuff = prevUnd.buffer(0, join_style=2)
                    intersectObj = selfBuff.intersection(prevUndBuff)
                    if intersectObj:
                        if intersectObj.geom_type == 'MultiPolygon':
                            for plygn in intersectObj:
                                intersectCoords.append(
                                        zip(plygn.exterior.xy[0],
                                        plygn.exterior.xy[1]))
                                prevSplitUndObj = (prevUndBuff.difference(
                                        intersectObj))
                                newSplitUndulationObj = selfBuff.difference(
                                        intersectObj)
                        elif intersectObj.geom_type == 'Polygon':
                            intersectCoords.append(
                                    zip(intersectObj.exterior.xy[0],
                                    intersectObj.exterior.xy[1]))
                            prevSplitUndObj = (prevUndBuff.difference(
                                    intersectObj))
                            newSplitUndulationObj = selfBuff.difference(
                                    intersectObj)

                        prevSplitUndulationCoords = self.splitObject(
                                prevSplitUndObj)
                        newSplitUndulationCoords = self.splitObject(
                                newSplitUndulationObj)
                        for coordinateSet3 in intersectCoords:
                            if self.clean(coordinateSet3):
                                print coordinateSet3
                                newUndulation = Undulation(
                                    coords=self.clean(coordinateSet3),
                                    layer=self.layer,
                                    angleLabel=(self.angle,
                                        prevUnd.angle), check=False)  
                            else: continue           
                        for coordinateSet4 in prevSplitUndulationCoords:
                            if self.clean(coordinateSet4):
                                prevSplitUndulation = Undulation(
                                    coords=self.clean(coordinateSet4),
                                    layer=self.layer,
                                    angleLabel=prevUnd.angle,
                                    check=False)             
                            else: continue
                        for coordinateSet5 in newSplitUndulationCoords:
                            if self.clean(coordinateSet5):
                                newSplitUndulation = Undulation(
                                    coords=self.clean(coordinateSet5), 
                                    layer=self.layer, 
                                    angleLabel=self.angle,
                                    check=False)            
                            else: continue
                        
                        Undulation._instances.remove(prevUnd)
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

# tape1 = Tape(angle=0)
# testResin2 = Resin(
#         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)])
# testResin2 = Resin(
#         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)])

# tape3 = Tape(angle=90)
# testResin3 = Resin(
#         coords=[(-6.0, 75.0), (-5.0, 75.0), (-5.0, -75.0), (-6.0, -75.0)])
# testResin4 = Resin(
#         coords=[(5.0, 75.0), (6.0, 75.0), (6.0, -75.0), (5.0, -75.0)])

# # tape4 = Tape(coords=[(6.0, 75.0), (16.0, 75.0), (16.0, -75.0), (6.0, -75.0)])
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


