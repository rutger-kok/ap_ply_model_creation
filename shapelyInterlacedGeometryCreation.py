from shapely.geometry import Polygon
from shapely import affinity
from shapely.geometry import Point
import numpy as np
from shapely.ops import cascaded_union
from itertools import count

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
            angleLabel=None):
        
        self.createObject(angle, coords, holes, layer, angleLabel)
        # first check intersection with own type e.g. for tapes check 
        # intersect with other tapes
        intersectSelfCoordList, splitSelfObjList = self.checkIntersect(
                list(self.getinstances(self.layer)))
        if intersectSelfCoordList:
            for intersectSelfCoords in intersectSelfCoordList:
                intersectSelfObj = self.newInstance(
                        tempCoords=intersectSelfCoords, tempLayer=self.layer+1,
                        tempAngleLabel=self.angle)
        
        intersectOtherCoordList, splitOtherObjList = self.checkIntersect(
                list(Resin.getinstances(self.layer)))
        if intersectOtherCoordList:
            for intersectOtherCoords in intersectOtherCoordList:
                intersectOtherObj = Resin(
                        coords=intersectOtherCoords, layer=self.layer+1,
                        angleLabel=self.angle)
                undulationObj = Undulation(
                        coords=intersectOtherCoords, layer=self.layer,
                        angleLabel=self.angle)
        if splitOtherObjList and splitSelfObjList:
            totalSplitObjList = splitSelfObjList.intersection(splitOtherObjList)
            for splitCoords in self.splitObject(totalSplitObjList):
                splitObj = self.newInstance(
                        tempCoords=splitCoords, tempLayer=self.layer,
                        tempAngleLabel=self.angle)
        elif splitOtherObjList:
            for splitCoords in self.splitObject(splitOtherObjList):
                splitObj = self.newInstance(
                        tempCoords=splitCoords, tempLayer=self.layer,
                        tempAngleLabel=self.angle)
        elif splitSelfObjList:
            for splitCoords in self.splitObject(splitSelfObjList):
                splitObj = self.newInstance(
                        tempCoords=splitCoords, tempLayer=self.layer,
                        tempAngleLabel=self.angle)

        if not intersectSelfCoordList and not intersectOtherCoordList:
            for pgon in self._instances:
                if self.almost_equals(pgon, decimal=3):
                    break
            else:
                self._instances.append(self)

    def createObject(self, angle=0, coords=None , holes=None, layer=1,
            angleLabel=None):
        if coords == None:
            coords = (
                    [(-100.0, -w/2), (100.0, -w/2), (100.0, w/2), (-100.0, w/2)]
                    )
        
        # create initial polygon, rotate it, then initialize the tape object
        a = Polygon(coords, holes)
        b = affinity.rotate(a, angle, origin=centerpoint)
        c = b.intersection(boundary) # trim tape using boundary
        new_coords = zip(c.exterior.xy[0], c.exterior.xy[1])
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
                mergedObj = cascaded_union(objectList).buffer(-1*10**-8)
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
                        intersectObjBuffered = intersectObj.buffer(1*10**-10)
                        differenceObj = self.difference(intersectObjBuffered)
                elif intersectObj.geom_type == 'Polygon':
                    intersectCoords.append(
                            zip(intersectObj.exterior.xy[0], 
                            intersectObj.exterior.xy[1]))
                    intersectObjBuffered = intersectObj.buffer(1*10**-10)
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
        return 'Layer: {}, Number: {}'.format(self.layer, self.objNumber)
    def __str__(self):
        return 'Layer: {}, Number: {}'.format(self.layer, self.objNumber)               

# the Resin class is a subclass of the Tape class. Aside from object creation, 
# the class defines its own __init__ method as its intersection behaviour is
# different (resin intersecting tape -> no undulation). 
class Resin(Tape):
    _ids = count(0)
    _instances = []

    def __init__(self, angle=0, coords=None , holes=None, layer=1,
            angleLabel=None):
        self.createObject(angle, coords, holes, layer, angleLabel)
        # first check intersection with own type e.g. for tapes check intersect
        # with other tapes
        intersectSelfCoordList, splitSelfObjList = self.checkIntersect(
                list(self.getinstances(self.layer)))
        if intersectSelfCoordList:
            for intersectSelfCoords in intersectSelfCoordList:
                intersectSelfObj = self.newInstance(
                        tempCoords=intersectSelfCoords, tempLayer=self.layer+1,
                        tempAngleLabel=self.angle)

        # check intersect with other e.g. for tapes check intersect with resin 
        intersectOtherCoordList, splitOtherObjList = self.checkIntersect(
                list(Tape.getinstances(self.layer)))
        if intersectOtherCoordList:
            for intersectOtherCoords in intersectOtherCoordList:
                intersectOtherObj = self.newInstance(
                        tempCoords=intersectOtherCoords,
                        tempLayer=self.layer+1, tempAngleLabel=self.angle)

        if splitOtherObjList and splitSelfObjList:
            totalSplitObjList = splitSelfObjList.intersection(splitOtherObjList)
            for splitCoords in self.splitObject(totalSplitObjList):
                splitObj = self.newInstance(
                        tempCoords=splitCoords, tempLayer=self.layer,
                        tempAngleLabel=self.angle)
        elif splitOtherObjList:
            for splitCoords in self.splitObject(splitOtherObjList):
                splitObj = self.newInstance(
                        tempCoords=splitCoords, tempLayer=self.layer,
                        tempAngleLabel=self.angle)
        elif splitSelfObjList:
            for splitCoords in self.splitObject(splitSelfObjList):
                splitObj = self.newInstance(
                        tempCoords=splitCoords, tempLayer=self.layer,
                        tempAngleLabel=self.angle)

        if not intersectSelfCoordList and not intersectOtherCoordList:
            for pgon in self._instances:
                if self.almost_equals(pgon, decimal=3):
                    break
            else:
                self._instances.append(self)

# the Undulation class is a subclass of the Tape class. The Undulation class
# defines its own __init__ method which defines how the undulation regions
# are recorded.
class Undulation(Tape):
    _ids = count(0)
    _instances = []
    def __init__(self, angle=0, coords=None , holes=None, layer=1,
            angleLabel=None, check=True):
        self.createObject(angle, coords, holes, layer, angleLabel)
        splitUndulationCoordsList = None
        if check:
            intersectSelfCoordList, splitSelfObjList = self.checkIntersect(
                    list(self.getinstances(self.layer)), buffer=False)
            if intersectSelfCoordList:
                for intersectSelfCoords in intersectSelfCoordList:
                    intersectSelfObj = self.newInstance(
                            tempCoords=intersectSelfCoords,
                            tempLayer=self.layer+1, tempAngleLabel=self.angle)
                    splitUndulationCoordsList = self.splitObject(
                            splitSelfObjList)

            else:
                if self.layer == 1:
                    newUndulation = Undulation(
                            coords=self.coordinates, layer=2,
                            angleLabel=self.angle, check=False)
                    bottomLayer = 1
                else:
                    newUndulation = Undulation(
                            coords=self.coordinates, layer=self.layer,
                            angleLabel=self.angle, check=False)
                    bottomLayer = self.layer - 1

                for prevUndulation in list(self.getinstances(bottomLayer)):
                    if self.intersection(prevUndulation):
                        prevSplitUndulationObj = None
                        newSplitUndulationObj = None
                        intersectCoords = []
                        intersectObj = self.intersection(prevUndulation)
                        if intersectObj:
                            if intersectObj.geom_type == 'MultiPolygon':
                                for plygn in intersectObj:
                                    intersectCoords.append(
                                            zip(plygn.exterior.xy[0],
                                            plygn.exterior.xy[1]))
                                    intersectObjBuffered = intersectObj.buffer(
                                            1*10**-10)
                                    prevSplitUndulationObj = (
                                            prevUndulation.difference(
                                                intersectObjBuffered))
                                    newSplitUndulationObj = self.difference(
                                            intersectObjBuffered)
                            elif intersectObj.geom_type == 'Polygon':
                                intersectCoords.append(
                                        zip(intersectObj.exterior.xy[0],
                                        intersectObj.exterior.xy[1]))
                                intersectObjBuffered = intersectObj.buffer(
                                        1*10**-10)
                                prevSplitUndulationObj = (
                                        prevUndulation.difference(
                                            intersectObjBuffered))
                                newSplitUndulationObj = self.difference(
                                        intersectObjBuffered)

                            prevSplitUndulationCoords = self.splitObject(
                                    prevSplitUndulationObj)
                            newSplitUndulationCoords = self.splitObject(
                                    newSplitUndulationObj)
                            for coordinateSet in intersectCoords:
                                newUndulation = Undulation(
                                        coords=coordinateSet, layer=bottomLayer,
                                        angleLabel=(self.angle,
                                            prevUndulation.angle), check=False)             
                            for coordinateSet2 in prevSplitUndulationCoords:
                                prevSplitUndulation = Undulation(
                                        coords=coordinateSet2,
                                        layer=bottomLayer,
                                        angleLabel=prevUndulation.angle,
                                        check=False)             
                            for coordinateSet3 in newSplitUndulationCoords:
                                newSplitUndulation = Undulation(
                                        coords=coordinateSet3, 
                                        layer=bottomLayer, 
                                        angleLabel=self.angle,
                                        check=False)             
                            
                            self._instances.remove(prevUndulation)
                    else:
                        continue
                else:
                    newUndulation = Undulation(
                            coords=self.coordinates, layer=bottomLayer,
                            angleLabel=self.angle, check=False)
        
        if splitUndulationCoordsList:
            for splitUndulationCoords in splitUndulationCoordsList:
                splitUndulationBottom = Undulation(
                        coords=splitUndulationCoords, layer=self.layer,
                        angleLabel=self.angle, check=False)     
                splitUndulationTop = Undulation(
                        coords=splitUndulationCoords, layer=self.layer+1,
                        angleLabel=self.angle, check=False)

        if not check:
            self._instances.append(self)


testResin2 = Resin(
        coords=[(-100.0, -6), (100.0, -6), (100.0, -5), (-100.0, -5)])
testResin2 = Resin(
        coords=[(-100.0, 6), (100.0, 6), (100.0, 5), (-100.0, 5)])
testResin3 = Resin(
        coords=[(-6.0, 75.0), (-5.0, 75.0), (-5.0, -75.0), (-6.0, -75.0)])
testResin4 = Resin(
        coords=[(5.0, 75.0), (6.0, 75.0), (6.0, -75.0), (5.0, -75.0)])


tape1 = Tape(angle=0)
tape3 = Tape(angle=90)
tape2 = Tape(angle=45)

# -----------------------------------------------------------------------------
# Plotting of regions
# import matplotlib.pyplot as plt

# f, axes = plt.subplots(3, 3, sharex='col', sharey='row')
# f.suptitle('Geometry per Layer')

# g = 1
# for i in range(3):
#     for j in range(3):
#         if list(Undulation.getinstances(g)):
#             for k in range(len(list(Undulation.getinstances(g)))):
#                 xi,yi = list(Undulation.getinstances(g))[k].exterior.xy
#                 axes[i,j].plot(xi,yi)
#                 axes[i,j].set_title('Layer: {}'.format(g))
#             xb, yb = boundary.exterior.xy
#             axes[i,j].plot(xb,yb)
#         else:
#             break
#         g += 1

# g = 1
# for i in range(3):
#     for j in range(3):
#         if list(Tape.getinstances(g)):
#             for k in range(len(list(Tape.getinstances(g)))):
#                 xi,yi = list(Tape.getinstances(g))[k].exterior.xy
#                 axes[i,j].plot(xi,yi)
#                 axes[i,j].set_title('Layer: {}'.format(g))
#             xb, yb = boundary.exterior.xy
#             axes[i,j].plot(xb,yb)
#         else:
#             break
#         g += 1
# plt.show(f)

