from sys import path
path.append(r'C:\Python27\Lib\site-packages')
from abaqus import *
from abaqusConstants import *
import mesh
import regionToolset
import math
import job
import displayGroupMdbToolset as dgm
from shapely.geometry import *
from shapely.geometry.polygon import Polygon

from shapely.geometry import Polygon
from shapely import affinity
from shapely.geometry import Point
import numpy as np
from shapely.ops import cascaded_union
from itertools import count
from polygonv7 import *

import os

session.viewports['Viewport: 1'].setValues(displayedObject=None)

mdb.models.changeKey(fromName='Model-1', toName='Layer 1')
laminateModel = mdb.models['Layer 1']
laminateAssembly = laminateModel.rootAssembly

t = 2.0  # layer thickness 

def definePart(objectList):
    for num, obj in enumerate(objectList):
        x0, y0 = obj.exterior.xy
        x0r, y0r = map(lambda x: round(x, 5), x0), map(lambda x: round(x, 5), y0)
        abqCoords = zip(x0r, y0r)

        profileSketch = mdb.models['Layer 1'].ConstrainedSketch(name='Object {}'.format(num), 
                sheetSize=200.0)
        for ind in range(len(abqCoords)-1):
            profileSketch.Line(point1=(abqCoords[ind][0], abqCoords[ind][1]), point2=(abqCoords[ind+1][0], abqCoords[ind+1][1]))
        
        objPart = laminateModel.Part(name='Object {}'.format(num), dimensionality=THREE_D, 
                type=DEFORMABLE_BODY)
        objPart.BaseSolidExtrude(sketch=profileSketch, depth=t)


test = definePart(list(Tape.getinstances(1)))