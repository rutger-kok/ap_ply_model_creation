from odbAccess import *
import numpy as np

path =r"C:\Workspace\impact\composite\DWT_XP\DWT_XP.odb"
odb = openOdb(path)

numFrames = len(odb.steps['Loading Step'].frames)
numElems = len(odb.rootAssembly.instances['Specimen Instance'].elementSets['Layer 1'].elements)
nodeLabels = {node.label: i for i,node in enumerate(odb.rootAssembly.instances['Specimen Instance'].nodes)}
elemLabels = {elem.label: j for j,elem in enumerate(odb.rootAssembly.instances['Specimen Instance'].elements)}
bottomLayerCells = odb.rootAssembly.instances['Specimen Instance'].elementSets['Layer 1'].elements

csv_array = np.zeros(shape=(numFrames*numElems,6))
counter = 0
for frame in odb.steps['Loading Step'].frames:
    strain = frame.fieldOutputs['LE']
    timeStep = frame.frameValue
    for cell in bottomLayerCells:
        elemNum = cell.label
        elemIndex = elemLabels[elemNum]
        strain11 = strain.values[elemIndex].data[0]
        strain22 = strain.values[elemIndex].data[1]
        nodeCoords = np.zeros(shape=(len(cell.connectivity),3))
        for m,nodeNum in enumerate(cell.connectivity):
            index = nodeLabels[nodeNum]
            x,y,z = strain.values[0].instance.nodes[index].coordinates
            nodeCoords[m] = x,y,z
        centroid = np.mean(nodeCoords, axis=0)
        csv_array[counter] = timeStep,centroid[0],centroid[1],centroid[2],strain11,strain22
        counter += 1

# -------------------------------------------------------------------
# Write to csv
csvPath = 'C:\\Workspace\\impact\\composite\\DWT_XP\\DWT_XP_strainContours.csv'
heading = 'time,x,y,z,e11,e22'
np.savetxt(csvPath,csv_array,delimiter=',')

# -------------------------------------------------------------------
# Extract force, displacement, velocity data

impactorRF = odb.rootAssembly.instances['Impactor Instance'].nodeSets['Impactor Reference Point']

csv2_array = np.zeros(shape=(numFrames,5))
counter2 = 0
for frame in odb.steps['Loading Step'].frames:
    vel = frame.fieldOutputs['V']
    disp = frame.fieldOutputs['U']
    acc = frame.fieldOutputs['A']
    timeStep = frame.frameValue

    subset_vel = vel.getSubset(region=impactorRF)
    subset_disp = disp.getSubset(region=impactorRF)
    subset_acc = acc.getSubset(region=impactorRF)

    vel_data = subset_vel.values[0].data[1] # velocity in 2 direction
    disp_data = subset_disp.values[0].data[1] # velocity in 2 direction
    acc_data = subset_acc.values[0].data[1] # velocity in 2 direction
    force_data = acc_data*5.0 # impactor mass = 5.0 kg

    csv2_array[counter2] = timeStep,disp_data,vel_data,acc_data,force_data
    counter2 += 1

# -------------------------------------------------------------------
# Write to csv
csv2Path = 'C:\\Workspace\\impact\\composite\\DWT_XP\\DWT_XP_forceData.csv'
np.savetxt(csv2Path,csv2_array,delimiter=',')



# data = zip(e11_avg,s11_avg)
# plotName = 'S11 vs LE11'
# xAxisTitle = 'LE11'
# yAxisTitle = 'S11'
# xyData = session.XYData(plotName, data)
# curve = session.Curve(xyData)
# xyPlot = session.XYPlot(plotName)
# chart = xyPlot.charts.values()[0]
# chart.setValues(curvesToPlot=(plotName, ))
# chart.axes1[0].axisData.setValues(useSystemTitle=False,
#                                     title=xAxisTitle)
# chart.axes2[0].axisData.setValues(useSystemTitle=False,
#                                     title=yAxisTitle)
# # Display the XY Plot in the current viewport
# vp = session.viewports[session.currentViewportName]
# vp.setValues(displayedObject=xyPlot)
    