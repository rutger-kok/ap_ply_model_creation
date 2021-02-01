# cut outs from scripts which may prove to be useful

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

# def surfacePairs():
# allInstances = laminateAssembly.instances.values()
# laminateAssembly.InstanceFromBooleanMerge(name='Merge All',
#         instances=allInstances, keepIntersections=ON, originalInstances=DELETE, 
#         mergeNodes=ALL, nodeMergingTolerance=1e-06, domain=BOTH)
#     allFaces = laminateAssembly.instances['Merge All-1'].faces[:]
#     pairs = {}
#     for face in allFaces:
#         currentFaceID = face.index
#         currentCoords = face.pointOn[0]
#         adjFaceObjs = face.getAdjacentFaces()
#         adjFaceIDs = [f.index for f in adjFaceObjs]
#         adjFaceCoords = [f.pointOn[0] for f in adjFaceObjs]
#         adjFaces = zip(adjFaceIDs,adjFaceCoords)
#         for (adjFID, adjFCoords) in adjFaces:
#             if (currentFaceID, adjFID) not in pairs.keys() and (adjFID, currentFaceID) not in pairs.keys():
#                 pairs[(currentFaceID,adjFID)] = (currentCoords, adjFCoords)
#             else: continue 
#     return pairs

# def createCohesivePairs():
# for (pairIDs, pairCoords) in pp.iteritems():
#     face1 = [inst.faces.findAt((pairCoords[0], )) for inst
#             in laminateAssembly.instances.values() 
#             if inst.faces.findAt(pairCoords[0], )]
#     face2 = [inst.faces.findAt((pairCoords[1], )) for inst
#             in laminateAssembly.instances.values() 
#             if inst.faces.findAt(pairCoords[1], )]
#     faceSet = []
#     for ff in face1+face2:
#         if ff not in faceSet:
#             faceSet.append(ff)
#     faceRegion1 = laminateAssembly.Surface(side1Faces=faceSet[0], name='Face {}'.format(pairIDs[0]))
#     faceRegion2 = laminateAssembly.Surface(side1Faces=faceSet[1], name='Face {}'.format(pairIDs[1]))
#     if len(faceSet) > 2:
#         print 'Error - more than 2 faces in pair'
#     laminateModel.SurfaceToSurfaceContactExp(name ='Coh {}-{}'.format(pairIDs[0], pairIDs[1]), 
#         createStepName='Initial', master = faceRegion1, slave = faceRegion2, 
#         mechanicalConstraint=KINEMATIC, sliding=FINITE, 
#         interactionProperty='Cohesive Surface', initialClearance=OMIT, 
#         datumAxis=None, clearanceRegion=None)
