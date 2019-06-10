from shapely.geometry import Polygon
from shapely import affinity
from shapely.geometry import Point
import numpy as np
from shapely.ops import cascaded_union
from itertools import count
from itertools import combinations
from shapely.geometry import shape, JOIN_STYLE
# define outside boundary of the impact specimen
boundary = Polygon([(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])


# import matplotlib.pyplot as plt
# test1 = Polygon([(7.0, -75.0), (17.0, -75.0), (17.0, 75.0), (7.0, 75.0)])
# xx, yy = test1.exterior.xy
# test2 = Polygon([(-6.00000002, -6.00000001), (-13.07106783, -6.00000001), (-6.00000002, 1.0710678), (-6.00000002, -6.00000001)])
# xx2, yy2 = test2.exterior.xy
# test3 = Polygon([(1.0710678, -6.00000002), (-6.00000001, -13.07106783), (-6.00000001, -6.00000002), (1.0710678, -6.00000002)])
# xx3, yy3 = test3.exterior.xy
# test4 = Polygon([(5.0, 4.99999999), (5.0, 4.99999995), (5.00000001, 4.99999995), (5.00000001, 4.99999992), (5.00000001, 4.99999992), (5.00000001, 4.99999989), (5.00000001, 4.99999989), (5.00000001, -2.07106782), (2.07106784, -4.99999999), (-5.0, -4.99999999), (-5.0, -4.99999995), (-5.00000001, -4.99999995), (-5.00000001, -4.99999992), (-5.00000001, -4.99999992), (-5.00000001, -4.99999989), (-5.00000001, -4.99999989), (-5.00000001, 2.07106782), (-2.07106784, 4.99999999), (5.0, 4.99999999)])
# xx4, yy4 = test4.exterior.xy
# test5 = Polygon([(-1.0710678, 6.00000002), (6.00000001, 13.07106783), (6.00000001, 6.00000002), (-1.0710678, 6.00000002)])
# xx5, yy5 = test5.exterior.xy
# test6 = Polygon([(13.07106783, 6.00000001), (6.00000002, -1.0710678), (6.00000002, 6.00000001), (13.07106783, 6.00000001)])
# xx6, yy6 = test6.exterior.xy
# test7 = Polygon([(50.0, 42.92893219), (13.07106782, 6.00000001), (6.00000001, 6.00000001), (6.00000001, 13.07106782), (50.0, 57.07106781), (50.0, 42.92893219)])
# xx7, yy7 = test7.exterior.xy
# xb, yb = boundary.exterior.xy
# plt.plot(xx,yy)
# plt.plot(xx2,yy2)
# plt.plot(xx3,yy3)
# plt.plot(xx4,yy4)
# plt.plot(xx5,yy5)
# plt.plot(xx6,yy6)
# plt.plot(xx7,yy7)
# plt.plot(xb,yb)
# plt.show()

# define tape width and thickness
w = 10.0
t = 0.2

def splitObject(differenceObject):
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


def clean(coordinates):
    poly = Polygon(coordinates).buffer(0)
    tol = 1*10**-8# distance
    cf = 1.5  # cofactor
    cleaned = poly.buffer(-tol,join_style=2).buffer(
            tol*cf, join_style=2).intersection(poly).simplify(1*10**-10)
    try:
        new_coords = list(
                zip(cleaned.exterior.xy[0],cleaned.exterior.xy[1]))
        return new_coords
    except AttributeError:
        return None
        

class machinePass():
    def __init__(self, angle=0, coords=None, resinW=1):

        self.tapeInstances = []
        self.otherInstances = []

        if coords == None:
            coords = (
                    [(-100.0, -w/2), (100.0, -w/2), (100.0, w/2), (-100.0, w/2)]
                    )
        # create Tape
        Tape(self, coords=coords, angle=angle)

        tapePath = sorted(
                self.tapeInstances,
                key = lambda x: (round(x.centroid.coords[0][0],6),
                    round(x.centroid.coords[0][1],6)))

        for index in range(len(tapePath)-1):
            tapeInstance1 = tapePath[index]
            tapeInstance2 = tapePath[index+1]
            if tapeInstance1.layer != tapeInstance2.layer:
                for plygn in self.otherInstances:
                    if plygn[0].buffer(1*10**-8).intersects(tapeInstance1) and plygn[0].buffer(1*10**-8).intersects(tapeInstance2):
                        x0, y0 = plygn[0].exterior.xy
                        coordinateSet = zip(x0,y0)
                        if clean(coordinateSet):
                            newUndBottomLayer = Undulation(
                                    coords= clean(coordinateSet),
                                    layer=plygn[1]-1, angleLabel=angle)   
                            newUndTopLayer = Undulation(
                                    coords=clean(coordinateSet),
                                    layer=plygn[1], angleLabel=angle)
                        else: continue 

                        validResinObject = map(lambda x: x.buffer(0), Resin.getinstances(plygn[1]-1))
                        #iterate over resin objects
                        for resinObject in validResinObject:
                            diffObject = plygn[0].buffer(
                                    0, join_style=2)
                            if diffObject.intersection(resinObject):
                                Resin._instances.remove(resinObject)
                                resinSplitObjectList = resinObject.difference(
                                        diffObject)
                                resinSplitCoords = splitObject(
                                        resinSplitObjectList)       

                                # create split resin regions
                                for coordinateSet2 in resinSplitCoords:
                                    if clean(coordinateSet2):
                                        splitResin = Resin(
                                            coords=coordinateSet2,
                                            layer=plygn[1]-1,
                                            angleLabel=angle, check=False)
                                    else: continue

                    else: continue
            else:
                for plygn in self.otherInstances:
                    if plygn[0].buffer(1*10**-8).intersects(tapeInstance1) and plygn[0].buffer(1*10**-8).intersects(tapeInstance2):
                        newTapeObj = cascaded_union([tapeInstance1.buffer(1*10**-8),tapeInstance2.buffer(1*10**-8),plygn[0].buffer(1*10**-8)])
                        x1,y1 = newTapeObj.exterior.xy
                        newTapeCoords = zip(x1,y1)
                        if clean(newTapeCoords):
                            newTape = Tape(self, coords=clean(newTapeCoords), layer=tapeInstance1.layer, check=False)
                            Tape._instances.remove(tapeInstance1)
                            Tape._instances.remove(tapeInstance2)
                        else: 
                            print 'cant clean'
        

        holla=[(coords[0][0], coords[0][1]),
                     (coords[1][0], coords[1][1]),
                     (coords[1][0],coords[1][1]-resinW),
                     (coords[0][0], coords[0][1]-resinW)]
        print holla
        yo = [(coords[2][0], coords[2][1]),
                     (coords[3][0], coords[3][1]),
                     (coords[3][0],coords[3][1]+resinW),
                     (coords[2][0], coords[2][1]+resinW)]
        print yo
        
        # create accompanying resin regions
        resinCoordsTop = Resin(
                coords=[(coords[0][0], coords[0][1]),
                     (coords[1][0], coords[1][1]),
                     (coords[1][0],coords[1][1]-resinW),
                     (coords[0][0], coords[0][1]-resinW)], angle=angle)

        resinCoordsBottom = Resin(
                coords=[(coords[2][0], coords[2][1]),
                     (coords[3][0], coords[3][1]),
                     (coords[3][0],coords[3][1]+resinW),
                     (coords[2][0], coords[2][1]+resinW)], angle=angle)


# defines the Tape class, which is a subclass of the Shapely Polygon class.
class Tape(Polygon):
    _ids = count(0)  # counts the number of instances created
    _instances = []  # initializes a list of Tape instances 
    def __init__(self, passNum, angle=0, coords=None, layer=1,
            angleLabel=None, check=True, step1=True):

        self.passNumber = passNum
        self.createObject(angle, coords, layer, angleLabel)

        if check: 
            if step1:
                allObjInLayer = list(list(Tape.getinstances(self.layer)) + 
                    list(Undulation.getinstances(self.layer)) + 
                    list(Resin.getinstances(self.layer)))

                tapeIntAllCoordList, tapeSplitAll = self.checkIntersect(
                    allObjInLayer)
                if tapeIntAllCoordList:
                    for tapeIntAllCoords in tapeIntAllCoordList:
                        if clean(tapeIntAllCoords):
                            mergedTapeLayer = Tape(
                                passNum = self.passNumber,
                                coords=clean(tapeIntAllCoords),
                                layer=self.layer+1, angleLabel=self.angle)
                        else: continue
                    for tapeSplitAllCoords in splitObject(tapeSplitAll):
                        if clean(tapeSplitAllCoords):
                            tapeSplitAllObj = Tape(
                                passNum=self.passNumber,
                                coords=clean(tapeSplitAllCoords),
                                layer=self.layer, angleLabel=self.angle,
                                step1=False)     
                        else: continue           
                else:
                    if clean(self.coordinates):
                        Tape(coords=clean(self.coordinates),
                                passNum=self.passNumber,
                                layer=self.layer, angleLabel = self.angle,
                                step1=False)

            else:
                resinUndulationMerge = list(
                        list(Resin.getinstances(self.layer-1))+
                        list(Undulation.getinstances(self.layer-1)))
                tapeIntMergeCoordList, tapeSplitMerge = self.checkIntersect(
                        resinUndulationMerge)
                if tapeIntMergeCoordList:
                    passNum.otherInstances.extend([(Polygon(crds), self.layer) for crds in tapeIntMergeCoordList])

                # create split tape object (for non-intersecting Tape regions)
                if tapeSplitMerge:
                    for tapeSplitCoords in splitObject(tapeSplitMerge):
                        if clean(tapeSplitCoords):
                            tapeSplitObj = Tape(
                                passNum=self.passNumber,
                                coords=clean(tapeSplitCoords),
                                layer=self.layer,
                                angleLabel=self.angle, check=False)
                        else: continue

                # if the current Tape doesnt intersect with any other object it
                # is added to the instance list
                if not tapeIntMergeCoordList:
                    for pgon in list(self.getinstances(self.layer)):
                        if self.almost_equals(pgon, decimal=3):
                            break
                    else:
                        if self.area > 0.001:
                            self._instances.append(self)
                            passNum.tapeInstances.append(self)
        else:
            for pgon in list(self.getinstances(self.layer)):
                if self.almost_equals(pgon, decimal=3):
                    break
            else:
                if self.area > 0.001:
                    self._instances.append(self)
                    passNum.tapeInstances.append(self)

    def createObject(self, angle=0, coords=None ,layer=1,
            angleLabel=None, rotationPoint=None):

        a = Polygon(coords)
        b = affinity.rotate(a, angle, origin=Point(0,0))
        c = b.intersection(boundary) # trim tape using boundary
        x0r = map(lambda x: round(x, 8), c.exterior.xy[0]) 
        y0r = map(lambda x: round(x, 8), c.exterior.xy[1])
        new_coords = zip(x0r,y0r)
        # new_coords = zip(c.exterior.xy[0],c.exterior.xy[1])
        Polygon.__init__(self,new_coords) # initialize rotated tape

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

    # this method checks whether the current object intersects with any of the
    # objects in "objectList"
    def checkIntersect(self, objectList):
        differenceObj = None
        intersectCoords = []
        validObjectList = map(lambda x: x.buffer(0), objectList)
        if validObjectList:
            mergedObj = cascaded_union(validObjectList).buffer(1*10**-8, join_style=2)
            selfBuffer = self.buffer(1*10**-8, join_style=2)
            intersectObj = selfBuffer.intersection(mergedObj).buffer(0)
            if intersectObj and not intersectObj.is_empty:
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
                    differenceObj = self.buffer(0).difference(intersectObj)
                else:
                    print 'Not a polygon or multipolygon'     
        return intersectCoords, differenceObj

    def __repr__(self):
        return 'Tape: Layer: {}, Number: {}, Coordinates: {}'.format(
                self.layer, self.objNumber, self.coordinates)
    def __str__(self):
        return 'Tape: Layer: {}, Number: {}, Coordinates: {}'.format(
                self.layer, self.objNumber, self.coordinates)          

# the Resin class is a subclass of the Tape class. Aside from object creation, 
# the class defines its own __init__ method as its intersection behaviour is
# different (resin intersecting tape -> no undulation). 
class Resin(Tape):
    _ids = count(0)
    _instances = []

    def __init__(self, angle=0, coords=None , layer=1, angleLabel=None,
            check=True):

        self.createObject(angle, coords, layer, angleLabel)

        if check:
            allObjInLayer = list(list(Tape.getinstances(self.layer)) + 
                    list(Undulation.getinstances(self.layer)) + 
                    list(Resin.getinstances(self.layer)))

            resinIntAllCoordList, resinSplitAll = self.checkIntersect(
                    allObjInLayer)
            if resinIntAllCoordList:
                for resinIntAllCoords in resinIntAllCoordList:
                    if clean(resinIntAllCoords):
                        mergedResinLayer = Resin(
                            coords=clean(resinIntAllCoords),
                            layer=self.layer+1, angleLabel=self.angle)
                    else: continue
                for resinSplitAllCoords in splitObject(resinSplitAll):
                    if clean(resinSplitAllCoords):
                        resinSplitAllObj = Resin(
                            coords=clean(resinSplitAllCoords),
                            layer=self.layer, angleLabel=self.angle,
                            check=False)     
                    else: continue           
            else:
                if clean(self.coordinates):
                    Resin(coords=clean(self.coordinates), layer=self.layer, 
                        angleLabel = self.angle, check=False)

        else:
            if self.area > 0.001:
                self._instances.append(self)

    def __repr__(self):
        return 'Resin: Layer: {}, Number: {}, Coordinates: {}'.format(
                self.layer, self.objNumber, self.coordinates)
    def __str__(self):
        return 'Resin: Layer: {}, Number: {}, Coordinates: {}'.format(
                self.layer, self.objNumber, self.coordinates)      

# the Undulation class is a subclass of the Tape class. The Undulation class
# defines its own __init__ method which defines how the undulation regions
# are recorded.
class Undulation(Tape):
    _ids = count(0)
    _instances = []
    def __init__(self, angle=0, coords=None, layer=1,
            angleLabel=None, check=True):
        self.createObject(angle, coords, layer, angleLabel)
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

                        prevSplitUndulationCoords = splitObject(
                                prevSplitUndObj)
                        newSplitUndulationCoords = splitObject(
                                newSplitUndulationObj)
                        for coordinateSet3 in intersectCoords:
                            if clean(coordinateSet3):
                                newUndulation = Undulation(
                                    coords=clean(coordinateSet3),
                                    layer=self.layer,
                                    angleLabel=(self.angle,
                                        prevUnd.angle), check=False)  
                            else: continue           
                        for coordinateSet4 in prevSplitUndulationCoords:
                            if clean(coordinateSet4):
                                prevSplitUndulation = Undulation(
                                    coords=clean(coordinateSet4),
                                    layer=self.layer,
                                    angleLabel=prevUnd.angle,
                                    check=False)             
                            else: continue
                        for coordinateSet5 in newSplitUndulationCoords:
                            if clean(coordinateSet5):
                                newSplitUndulation = Undulation(
                                    coords=clean(coordinateSet5), 
                                    layer=self.layer, 
                                    angleLabel=self.angle,
                                    check=False)            
                            else: continue
                        
                        Undulation._instances.remove(prevUnd)
                else:
                    continue
            else:
                if clean(self.coordinates):
                    newUndulation = Undulation(
                        coords=clean(self.coordinates), layer=self.layer,
                        angleLabel=self.angle, check=False)

        else:
            for pgon in list(self.getinstances(self.layer)):
                if self.almost_equals(pgon, decimal=1):
                    continue
            else:
                if self.area > 0.001:
                    self._instances.append(self)

    def __repr__(self):
        return 'Undulation: Layer: {}, Number: {}, Coordinates: {}'.format(
                self.layer, self.objNumber, self.coordinates)
    def __str__(self):
        return 'Undulation: Layer: {}, Number: {}, Coordinates: {}'.format(
                self.layer, self.objNumber, self.coordinates)      



# -----------------------------------------------------------------------------
# This section of the script only runs if this script is run directly (i.e. as
# long as this script is not imported). It plots the geometries for testing
# purposes.

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # Plotting of regions

    # NOTE when specifying coordinates must move clockwise from top-left!

    # test case 1: two tapes and resin (touching) crossing tape and resin
    # PASS: mergebuffer: 1*10**-8, selfbuffer: 1*10**-8
    # testTape1 = machinePass(angle=0)
    # testTape2 = machinePass(angle=90)
    # testTape3 = machinePass(coords=[(-100.0, -7.0), (100.0, -7.0), (100.0, -17.0), (-100.0, -17.0)], angle=90)

    # test case 7: two angled tapes touching
    # PASS: mergebuffer: 1*10**-8, selfbuffer: 1*10**-8
    # testTape1 = Tape(coords=[(-100.0, -5.0), (100.0, -5.0), (100.0, 5.0), (-100.0, 5.0)], angle=45)
    # testTape2 = Tape(coords=[(-100.0, -5.0), (100.0, -5.0), (100.0, -15.0), (-100.0, -15.0)], angle=45, rotateP=testTape1.rotationPoint)

    # test case 8: angled tape crossing resin and tape
    # PASS: mergebuffer: 1*10**-8, selfbuffer: 1*10**-8
    # testTape1 = Tape(angle=0)
    # testResin1 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)])
    # testResin2 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)])
    # testTape2 = Tape(angle=45)

    # test case 9: angled tape and resin crossing resin and tape
    # PASS: mergebuffer: 1*10**-8, selfbuffer: 1*10**-8
    # testTape1 = Tape(angle=0)
    # testResin1 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)])
    # testResin2 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)])
    # testTape2 = Tape(angle=45)
    # testResin3 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=45, rotateP=testTape2.rotationPoint)
    # testResin4 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=45, rotateP=testTape2.rotationPoint)

    # test case 10: two angled tapes with touching resin regions crossing tape and resin region
    # PASS: mergebuffer: 1*10**-8, selfbuffer: 1*10**-8
    # testTape1 = Tape(angle=0)
    # testResin1 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)])
    # testResin2 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)])
    # testTape2 = Tape(angle=45)
    # testResin3 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=45, rotateP=testTape2.rotationPoint)
    # testResin4 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=45, rotateP=testTape2.rotationPoint)
    # testTape3 = Tape(coords=[(-100.0, -7.0), (100.0, -7.0), (100.0, -17.0), (-100.0, -17.0)], angle=45, rotateP=testTape2.rotationPoint)
    # testResin5 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -7.0), (-100.0, -7.0)], angle=45, rotateP=testTape2.rotationPoint)
    # testResin6 = Resin(
    #         coords=[(-100.0, -17.0), (100.0, -17.0), (100.0, -18.0), (-100.0, -18.0)], angle=45, rotateP=testTape2.rotationPoint)

    # test case 11: angled tape + resin regions intersecting angled tape + resin regions
    # PASS: mergebuffer: 1*10**-8, selfbuffer: 1*10**-8
    # testTape1 = Tape(angle=-45)
    # testResin1 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=-45, rotateP=testTape1.rotationPoint)
    # testResin2 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=-45, rotateP=testTape1.rotationPoint)
    # testTape2 = Tape(angle=45)
    # testResin3 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=45, rotateP=testTape2.rotationPoint)
    # testResin4 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=45, rotateP=testTape2.rotationPoint)

    # test case 12: two angled tape + (touching) resin regions intersecting angled tape + resin regions
    # PASS: mergebuffer: 1*10**-8, selfbuffer: 1*10**-8
    # testTape0 = machinePass(angle=0)
    # testTape1 = machinePass(angle=60)
    # testResin1 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=-45)
    # testResin2 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=-45)
    # testTape2 = machinePass(angle=120)
    # testResin3 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=45)
    # testResin4 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=45)
    # testTape3 = machinePass(coords=[(-100.0, -7.0), (100.0, -7.0), (100.0, -17.0), (-100.0, -17.0)], angle=45)
    # testResin5 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -7.0), (-100.0, -7.0)], angle=45)
    # testResin6 = Resin(
    #         coords=[(-100.0, -17.0), (100.0, -17.0), (100.0, -18.0), (-100.0, -18.0)], angle=45)

    # test case 13: two angled tapes + resin regions crossing two angled tapes + resin regions
    testTape1 = machinePass(angle=-45)
    # testResin1 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=-45)
    # testResin2 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=-45)
    testTape2 = machinePass(coords=[(-100.0, -17.0), (100.0, -17.0), (100.0, -7.0), (-100.0, -7.0)], angle=-45)
    # testResin3 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -7.0), (-100.0, -7.0)], angle=-45)
    # testResin4 = Resin(
    #         coords=[(-100.0, -17.0), (100.0, -17.0), (100.0, -18.0), (-100.0, -18.0)], angle=-45)
    testTape3 = machinePass(angle=45)
    # testResin5 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=45)
    # testResin6 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=45)
    testTape4 = machinePass(coords=[(-100.0, -17.0), (100.0, -17.0), (100.0, -7.0), (-100.0, -7.0)], angle=45)
    # testResin7 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -7.0), (-100.0, -7.0)], angle=45)
    # testResin8 = Resin(
    #         coords=[(-100.0, -17.0), (100.0, -17.0), (100.0, -18.0), (-100.0, -18.0)], angle=45)


    # testTape1 = machinePass(angle=0)
    # testResin1 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=0)
    # testResin2 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=0)
    # testTape2 = machinePass(angle=60)
    # testResin3 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=60)
    # testResin4 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=60)
    # testTape3 = machinePass(angle=120)
    # testResin5 = Resin(
    #         coords=[(-100.0, -6.0), (100.0, -6.0), (100.0, -5.0), (-100.0, -5.0)], angle=120)
    # testResin6 = Resin(
    #         coords=[(-100.0, 6.0), (100.0, 6.0), (100.0, 5.0), (-100.0, 5.0)], angle=120)

    # ------------------------------------------------------------------------

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
    
    print len(Tape._instances)
    print len(Resin._instances)
    print len(Undulation._instances)

    # print Tape._instances
    # print Resin._instances
    # print Undulation._instances

    # f4, axes = plt.subplots(3, 3, sharex='col', sharey='row')
    # f4.suptitle('Individual items per Layer')
    # g = 0
    # for i in range(3):
    #     for j in range(3):
    #         try:
    #             xi,yi = list(Resin.getinstances(1))[g].exterior.xy
    #             axes[i,j].plot(xi,yi)
    #             xb, yb = boundary.exterior.xy
    #             axes[i,j].plot(xb,yb)
    #             g += 1
    #         except IndexError:
    #             break

    # plt.show(f4)


