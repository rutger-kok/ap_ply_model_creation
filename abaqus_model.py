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

from abaqus import *
from abaqusConstants import *
from caeModules import *
import os


class Abaqus_Model():
    '''
    This class is used to create Abaqus models and contains a few useful
    methods for running and saving jobs, finding the faces of objects
    and so forth.
    '''
    def __init__(self, model_name, dir_name=None):
        self.model_name = model_name
        mdb.Model(name=self.model_name, modelType=STANDARD_EXPLICIT)
        self.model = mdb.models[self.model_name]
        self.assembly = self.model.rootAssembly
        # set the work directory (change this directory if necessary)
        if dir_name is None:
            wd = 'C:\\Workspace\\{}'.format(self.model_name)
        else:
            wd = 'C:\\Workspace\\{}\\{}'.format(dir_name, self.model_name)
        if not os.path.exists(wd):
            os.makedirs(wd)
        os.chdir(wd)

    def create_explicit_step(self, field_outputs=None, intervals=50,
                             step_name='Loading Step'):
        '''
        Method to create an Explicit step in an Abaqus simulation. This
        method also defines the desired field outputs and applies mass
        scaling to any elements within the model for which the stable
        time increment < 1e-6.
        '''
        # create step, 50 output intervals by default
        self.model.ExplicitDynamicsStep(
            name=step_name, previous='Initial', timePeriod=self.time,
            massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                          BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))
        # define default field outputs if they are not specified
        if field_outputs is None:
            field_outputs = ('S', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS',
                             'CSDMG', 'SDV', 'STATUS')
        self.model.FieldOutputRequest(
            name='F-Output-1', createStepName=step_name,
            variables=field_outputs, numIntervals=intervals)

    def create_implicit_step(self, field_outputs=None, initial_increment=0.1,
                             step_name='Loading Step'):
        '''
        Method to create an implicit step in an Abaqus simultion.
        '''
        self.model.StaticStep(
            name=step_name, previous='Initial', initialInc=initial_increment)
        # define default field outputs if not specified
        if field_outputs is None:
            field_outputs = ('S', 'E', 'U', 'RF', 'CSDMG', 'CSTRESS')
        self.model.FieldOutputRequest(
            name='F-Output-1', createStepName=step_name,
            variables=field_outputs)

    def get_faces_by_bounding_box(self, bounding_box, instances=None):
        '''
        This method searches instances in an assembly to determine which
        faces are contained within a bounding box. The standard Abaqus
        getByBoundingBox only operates on a single instance. This method
        finds faces within a bounding box for the instances specified.

        Args:
            bounding_box (list): A list of coodinates defining a bounding box
                in the form: (x_min, y_min, z_min, x_max, y_max, z_max)
            instances (list): list of the instances (not instance keys!) to be
                searched.

        Returns:
            face_seq (sequence of faces - Abaqus type): A sequence of the faces
                within the bounding box.
        '''
        if instances is None:
            search_instances = self.assembly.instances.values()
        else:
            search_instances = instances
        x1, y1, z1, x2, y2, z2 = bounding_box
        i = 0
        for inst in search_instances:
            if inst.faces.getByBoundingBox(x1, y1, z1, x2, y2, z2):
                faces = inst.faces.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                if i == 0:
                    face_seq = faces
                else:
                    face_seq += faces
                i += 1
        try:
            return face_seq
        except UnboundLocalError:
            print 'No faces within bounding box'
            return None

    def get_edges_by_bounding_box(self, bounding_box, instances=None):
        '''
        This method searches instances in an assembly to determine which
        edges are contained within a bounding box.

        Args:
            bounding_box (list): A list of coodinates defining a bounding box
                in the form: (x_min, y_min, z_min, x_max, y_max, z_max)
            instances (list): list of the instances (not instance keys!) to be
                searched.

        Returns:
            edge_seq (sequence of edges - Abaqus type): A sequence of the edges
                within the bounding box.
        '''
        if instances is None:
            search_instances = self.assembly.instances.values()
        else:
            search_instances = instances
        x1, y1, z1, x2, y2, z2 = bounding_box
        i = 0
        for inst in search_instances:
            if inst.edges.getByBoundingBox(x1, y1, z1, x2, y2, z2):
                edges = inst.edges.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                if i == 0:
                    edge_seq = edges
                else:
                    edge_seq += edges
                i += 1
        try:
            return edge_seq
        except UnboundLocalError:
            print 'No edges within bounding box'
            return None

    def get_cells_by_bounding_box(self, bounding_box, instances=None):
        '''
        This method searches instances in an assembly to determine which
        edges are contained within a bounding box.

        Args:
            bounding_box (list): A list of coodinates defining a bounding box
                in the form: (x_min, y_min, z_min, x_max, y_max, z_max)
            instances (list): list of the instances (not instance keys!) to be
                searched.

        Returns:
            cell_seq (sequence of cells - Abaqus type): A sequence of the cells
                within the bounding box.
        '''
        if instances is None:
            search_instances = self.assembly.instances.values()
        else:
            search_instances = instances
        x1, y1, z1, x2, y2, z2 = bounding_box
        i = 0
        for inst in search_instances:
            if inst.cells.getByBoundingBox(x1, y1, z1, x2, y2, z2):
                cells = inst.cells.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                if i == 0:
                    cell_seq = cells
                else:
                    cell_seq += cells
                i += 1
        try:
            return cell_seq
        except UnboundLocalError:
            print 'No cells within bounding box'
            return None

    def create_job(self, cpus=4, subroutine_path=None):
        '''
        Method to create an Abaqus job. By default the job definition
        specifies the use of 4 cpus in parallel.
        '''
        if subroutine_path:
            path = subroutine_path
        else:
            path = ''
        mdb.Job(
            name=self.model_name, model=self.model_name, description='',
            atTime=None, waitMinutes=0, waitHours=0, queue=None, type=ANALYSIS,
            memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
            nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
            contactPrint=OFF, historyPrint=OFF, userSubroutine=path,
            scratch='', resultsFormat=ODB, activateLoadBalancing=False,
            parallelizationMethodExplicit=DOMAIN, numDomains=cpus,
            multiprocessingMode=DEFAULT, numCpus=cpus, memory=90)

    def run_job(self):
        # submit job
        mdb.jobs[self.model_name].submit(consistencyChecking=OFF)
        # wait until job is complete
        mdb.jobs[self.model_name].waitForCompletion()

    def save_model(self):
        # save model
        mdb.saveAs(pathName=self.model_name)


if __name__ == '__main__':
    # Tests
    mdl = Abaqus_Model('TestModel')  # Create model
    mdl.time = 1.0  # set time attribute for testing
    mdl.create_explicit_step()
    mdl.create_job()
    mdl.save_model()
