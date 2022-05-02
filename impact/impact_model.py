'''
This module is part of a library used to generate AP-PLY composite
laminate geometries in Abaqus Explicit.
Copyright (C) 2022  Rutger Kok

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
USA
'''

import sys
# change these paths to point to your local Python installation package
# libraries and the AP-PLY model creation library.
sys.path.append('C:\\Python27\\Lib\\site-packages')
sys.path.append('C:\\Github\\ap_ply_model_creation')
from abaqus import *
from abaqusConstants import *
from shapely.geometry import Polygon
from ap_ply_3d import AP_PLY_3D
from ap_ply_elastic import AP_PLY_Elastic
import mesh
import regionToolset


class AP_PLY_Impact(AP_PLY_3D, AP_PLY_Elastic):
    def __init__(self, model_name, dir_name=None):
        AP_PLY_Elastic.__init__(self, model_name, dir_name)

    def set_test_parameters(self, time, output_intervals, energy,
                            shell_mesh, solid_mesh, symmetry_mode):
        '''Set drop weight tower test parameters'''
        self.time = time
        self.output_intervals = output_intervals
        self.energy = energy
        self.impactor_mass = 5.85  # kg
        # calculate initial impactor velocity
        self.init_velocity = -((2.0 * self.energy) / self.impactor_mass)**0.5
        # mesh densities
        self.shell_mesh = shell_mesh
        self.solid_mesh = solid_mesh
        self.symmetry_mode = symmetry_mode  # 'BC' or 'GEOM'

    def create_impactor_part(self):
        ''' Creates the hemispherical impactor for DWT simulations'''
        # Sketch the impactor (to create part by revolution)
        self.radius = 8.0
        sketch = self.model.ConstrainedSketch(
            name='Impactor Sketch', sheetSize=200.0)
        sketch.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
        sketch.ArcByCenterEnds(
            center=(0.0, 0.0), point1=(self.radius, 0.0),
            point2=(0.0, -self.radius), direction=CLOCKWISE)
        sketch.Line(point1=(0.0, -self.radius), point2=(0.0, 32.0))
        sketch.Line(point1=(0.0, 32.0), point2=(self.radius, 32.0))
        sketch.Line(point1=(self.radius, 32.0), point2=(self.radius, 0.0))
        self.impactor_part = self.model.Part(
            name='Impactor', dimensionality=THREE_D,
            type=DISCRETE_RIGID_SURFACE)
        self.impactor_part.BaseSolidRevolve(
            sketch=sketch, angle=360.0, flipRevolveDirection=OFF)
        rf_point_id = self.impactor_part.ReferencePoint(point=(0, 0, 0)).id
        rf_point = self.impactor_part.referencePoints[rf_point_id]
        rf_point_region = self.impactor_part.Set(
            referencePoints=(rf_point,), name='Impactor Reference Point')
        # Create shell
        cells = self.impactor_part.cells
        self.impactor_part.RemoveCells(cellList=cells[0:1])
        # Assign inertia/mass to impactor
        self.impactor_part.engineeringFeatures.PointMassInertia(
            name='Impactor Mass-Inertia', region=rf_point_region,
            mass=self.impactor_mass, alpha=0.0, composite=0.0)
        # Meshing
        elem_type_1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elem_type_2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        faces = self.impactor_part.faces
        picked_regions = (faces, )
        self.impactor_part.setElementType(
            regions=picked_regions, elemTypes=(elem_type_1, elem_type_2))
        self.impactor_part.seedPart(
            size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
        self.impactor_part.generateMesh()

    def create_clamp_part(self):
        '''Create the restraining clamps for DWT simulations'''
        coords = (44.0, 50.0)
        sketch = self.model.ConstrainedSketch(name='Clamp Sketch',
                                              sheetSize=200.0)
        sketch.CircleByCenterPerimeter(
            center=coords, point1=(coords[0] + 2.0, coords[1]))
        self.clamp_part = self.model.Part(
            name='Clamp', dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
        self.clamp_part.BaseSolidExtrude(sketch=sketch, depth=4.0)
        rf_point_id = self.clamp_part.ReferencePoint(
            point=(coords[0], coords[1], 0.0)).id
        rf_point = self.clamp_part.referencePoints[rf_point_id]
        rf_point_region = self.clamp_part.Set(
            referencePoints=(rf_point,), name='Clamp Reference Point')
        # Create shell
        cells = self.clamp_part.cells
        self.clamp_part.RemoveCells(cellList=cells[0:1])
        # Assign inertia/mass to clamps
        self.clamp_part.engineeringFeatures.PointMassInertia(
            name='Clamp Mass-Inertia', region=rf_point_region, mass=0.2,
            alpha=0.0, composite=0.0)
        # Mesh clamp
        faces = self.clamp_part.faces
        picked_regions = (faces, )
        elem_type_1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elem_type_2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        self.clamp_part.setElementType(
            regions=picked_regions, elemTypes=(elem_type_1, elem_type_2))
        self.clamp_part.seedPart(
            size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
        self.clamp_part.generateMesh()

    def create_bottom_plate(self):
        '''Create the plate on which specimens rest for DWT simulations'''
        sketch = self.model.ConstrainedSketch(name='Bottom Plate Sketch',
                                              sheetSize=200.0)
        sketch.rectangle(point1=(-75.0, -100.0), point2=(75.0, 100.0))
        sketch.rectangle(point1=(-37.5, -62.5), point2=(37.5, 62.5))
        self.b_plate_part = self.model.Part(
            name='Bottom Plate', dimensionality=THREE_D,
            type=DISCRETE_RIGID_SURFACE)
        self.b_plate_part.BaseShell(sketch=sketch)
        rf_point_id = self.b_plate_part.ReferencePoint(point=(0, 0, 0)).id
        rf_point = self.b_plate_part.referencePoints[rf_point_id]
        self.b_plate_part.Set(referencePoints=(rf_point,),
                              name='Bottom Plate Reference Point')
        # Do not need to assign intertia as plate is constrained in every DOF
        # Mesh bottom plate
        faces = self.b_plate_part.faces
        picked_regions = (faces, )
        elem_type_1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elem_type_2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        self.b_plate_part.setElementType(
            regions=picked_regions, elemTypes=(elem_type_1, elem_type_2))
        self.b_plate_part.seedPart(
            size=2.5, deviationFactor=0.1, minSizeFactor=0.1)
        self.b_plate_part.generateMesh()

    def create_assembly(self):
        # Create the part instances
        self.assembly.Instance(
            name='Impactor Instance', part=self.impactor_part, dependent=ON)
        self.assembly.Instance(
            name='Bottom Plate Instance', part=self.b_plate_part, dependent=ON)
        for n in range(1, 5):
            instance_name = 'Clamp Instance {}'.format(n)
            self.assembly.Instance(name=instance_name, part=self.clamp_part,
                                   dependent=ON)
        self.assembly.rotate(
            instanceList=('Bottom Plate Instance', 'Clamp Instance 1',
                          'Clamp Instance 2', 'Clamp Instance 3',
                          'Clamp Instance 4'),
            axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0),
            angle=-90.0)
        self.assembly.rotate(
            instanceList=('Bottom Plate Instance', 'Clamp Instance 1',
                          'Clamp Instance 2', 'Clamp Instance 3',
                          'Clamp Instance 4'),
            axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1.0, 0.0),
            angle=-90.0)
        self.assembly.translate(
            instanceList=('Impactor Instance', ),
            vector=(0.0, self.radius + self.l_thickness / 2.0 + 0.01, 0.0))
        self.assembly.translate(
            instanceList=('Clamp Instance 2', ), vector=(0.0, 0.0, -88.0))
        self.assembly.translate(
            instanceList=('Clamp Instance 3', ), vector=(-100.0, 0.0, 0.0))
        self.assembly.translate(
            instanceList=('Clamp Instance 4', ), vector=(-100.0, 0.0, -88.0))
        self.assembly.translate(
            instanceList=('Clamp Instance 1', 'Clamp Instance 2',
                          'Clamp Instance 3', 'Clamp Instance 4'),
            vector=(0.0, self.l_thickness / 2.0, 0.0))
        self.assembly.translate(
            instanceList=('Bottom Plate Instance', ),
            vector=(0.0, -self.l_thickness / 2.0, 0.0))

    def constrain_b_plate(self):
        # Encastre bottom plate
        instance = self.assembly.instances['Bottom Plate Instance']
        bc_region = instance.sets['Bottom Plate Reference Point']
        self.model.EncastreBC(name='Bottom Plate Encastre BC',
                              createStepName='Loading Step', localCsys=None,
                              region=bc_region)

    def constrain_clamps(self):
        # Fix clamps in all but y-direction
        all_instances = self.assembly.instances
        clamp_instances = [
            all_instances['Clamp Instance {}'.format(i)] for i in range(1, 5)]
        rf_points = [
            inst.referencePoints.values()[0] for inst in clamp_instances]
        clamp_region = self.assembly.Set(
            referencePoints=rf_points, name='Clamp BC Region')
        self.model.DisplacementBC(
            name='Clamp Displacement BC', createStepName='Loading Step',
            region=clamp_region, u1=0.0, u2=UNSET, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
        # Apply load to clamps
        self.model.SmoothStepAmplitude(
            name='Smoothing Amplitude', timeSpan=STEP,
            data=((0.0, 0.0), (1e-05, 1.0)))
        self.model.ConcentratedForce(
            name='Clamp Load BC', createStepName='Loading Step',
            region=clamp_region, cf2=-1.1, amplitude='Smoothing Amplitude',
            distributionType=UNIFORM, field='', localCsys=None)

    def constrain_impactor(self):
        # Fix impactor in all but y direction
        impactor_instance = self.assembly.instances['Impactor Instance']
        impactor_bc_region = impactor_instance.sets['Impactor Reference Point']
        self.model.DisplacementBC(
            name='Impactor Displacement BC', createStepName='Loading Step',
            region=impactor_bc_region, u1=0.0, u2=UNSET, u3=0.0, ur1=0.0,
            ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, fieldName='',
            distributionType=UNIFORM, localCsys=None)
        # Apply velocity BC to impactor
        self.model.Velocity(
            name='ImpactorVelocity', region=impactor_bc_region, field='',
            distributionType=MAGNITUDE, velocity1=0.0,
            velocity2=self.init_velocity, velocity3=0.0, omega=0.0)
        reg = impactor_instance.sets['Impactor Reference Point']
        self.model.HistoryOutputRequest(name='H-Output-2',
            createStepName='Loading Step', variables=('U2', 'V2', 'A2'),
            region=reg, sectionPoints=DEFAULT, rebar=EXCLUDE)


if __name__ == '__main__':

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    output_intervals = 50  # requested field output intervals
    energy = 30.0  # impact energy to simulate [J]
    symmetry_mode = 'GEOM'
    # Mesh sizes
    solid_mesh = 1.5
    shell_mesh = 1.5

    # Create model
    mdl = AP_PLY_Impact('DWT_Test')
    mdl.set_test_parameters(
        time, output_intervals, energy, shell_mesh, solid_mesh, symmetry_mode)
    mdl.l_thickness = 4.0  # for test purposes
    mdl.create_impactor_part()
    mdl.create_clamp_part()
    mdl.create_bottom_plate()
    mdl.create_explicit_step()
    mdl.create_assembly()
    mdl.constrain_b_plate()
    mdl.constrain_clamps()
    mdl.constrain_impactor()
    mdl.create_job()
    mdl.save_model()
