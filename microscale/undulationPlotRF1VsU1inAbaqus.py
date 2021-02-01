'''
Plot force vs displacement (in Abaqus) for a single undulation model.
'''

from odbAccess import *
import visualization
from caeModules import *

# Use the output database displayed in the current viewport
vp = session.viewports[session.currentViewportName]
odb = vp.displayedObject
if type(odb) != visualization.OdbType:
    raise ValueError, 'An odb must be displayed in the current viewport.'

# Alternatively use a path to an odb
# path = ''
# odb = openOdb(path)
# name = path[path.rfind('\\')+1:path.rfind('.')]

u1 = []
rf1 = []
for frame in odb.steps['Loading Step'].frames:
    force = frame.fieldOutputs['RF']
    rightCells = odb.rootAssembly.nodeSets['Right Coupling Region']
    subset_force = force.getSubset(region=rightCells)
    rf1_list = [subset_force.values[x].dataDouble for x in range(
                len(subset_force.values))]
    rf1.append(sum(rf1_list))
    displacement = frame.fieldOutputs['U']
    subset_displacement = displacement.getSubset(region=rightCells)
    u1_list = [subset_displacement.values[x].dataDouble[0] for x in range(
               len(subset_displacement.values))]
    u1.append((sum(u1_list) / len(u1_list)))

data = zip(u1, rf1)

# Plot data in Abaqus
fileName = odb.name[odb.name.rfind('/') + 1:odb.name.rfind('.')]
plotName = 'RF1 vs U1 - {}'.format(fileName)
xAxisTitle = 'U1'
yAxisTitle = 'RF1'
xyData = session.XYData(plotName, data)
curve = session.Curve(xyData)
xyPlot = session.XYPlot(plotName)
chart = xyPlot.charts.values()[0]
chart.setValues(curvesToPlot=(plotName, ))
chart.axes1[0].axisData.setValues(useSystemTitle=False,
                                  title=xAxisTitle)
chart.axes2[0].axisData.setValues(useSystemTitle=False,
                                  title=yAxisTitle)
# Display the XY Plot in the current viewport
vp = session.viewports[session.currentViewportName]
vp.setValues(displayedObject=xyPlot)
