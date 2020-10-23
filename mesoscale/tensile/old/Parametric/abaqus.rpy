# -*- coding: mbcs -*-
#
# Abaqus/CAE Release Unofficial Packaging Version replay file
# Internal Version: 2015_09_24-21.31.09 126547
# Run by s1342398 on Fri Apr 24 10:59:58 2020
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.55469, 1.55556), width=228.85, 
    height=154.311)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('interlacedRVEImplicitParametricMain.py', __main__.__dict__)
#: The model "0-90--20-5-0_18" has been created.
#: Model: C:/Workspace/3D_RVE/Parametric/0-90--20-5-0_18/0-90--20-5-0_18.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       1
#: Number of Node Sets:          120041
#: Number of Steps:              1
#: The model "0-90--25-1-0_18" has been created.
#: Model: C:/Workspace/3D_RVE/Parametric/0-90--25-1-0_18/0-90--25-1-0_18.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       1
#: Number of Node Sets:          22041
#: Number of Steps:              1
#: The model "0-90--25-2-0_18" has been created.
#: Model: C:/Workspace/3D_RVE/Parametric/0-90--25-2-0_18/0-90--25-2-0_18.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       1
#: Number of Node Sets:          48041
#: Number of Steps:              1
#: The model "0-90--25-3-0_18" has been created.
#: Model: C:/Workspace/3D_RVE/Parametric/0-90--25-3-0_18/0-90--25-3-0_18.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       1
#: Number of Node Sets:          84041
#: Number of Steps:              1
#: The model "0-90--25-4-0_18" has been created.
#: Model: C:/Workspace/3D_RVE/Parametric/0-90--25-4-0_18/0-90--25-4-0_18.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       1
#: Number of Node Sets:          130041
#: Number of Steps:              1
#: The model "0-90--25-5-0_18" has been created.
#: Model: C:/Workspace/3D_RVE/Parametric/0-90--25-5-0_18/0-90--25-5-0_18.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       1
#: Number of Node Sets:          186041
#: Number of Steps:              1
print 'RT script done'
#: RT script done
