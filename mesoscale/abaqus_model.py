'''
Module to create a Standard/Explicit Abaqus model

(c) Rutger Kok, 24/11/2020
'''

from abaqus import *
from abaqusConstants import *
from caeModules import *
import os


class AbaqusModel():
    def __init__(self, model_name, dir_name=None):
        self.model_name = model_name
        mdb.Model(name=self.model_name, modelType=STANDARD_EXPLICIT)
        self.model = mdb.models[self.model_name]
        self.assembly = self.model.rootAssembly
        # set the work directory
        if dir_name is None:
            wd = 'C:\\Workspace\\{}'.format(self.model_name)
        else:
            wd = 'C:\\Workspace\\{}\\{}'.format(dir_name, self.model_name)
        if not os.path.exists(wd):
            os.makedirs(wd)
        os.chdir(wd)

    def create_explicit_step(self, field_outputs=None, intervals=50,
                             step_name='Loading Step'):
        # mass scaling applied if increment < 1e-6
        self.model.ExplicitDynamicsStep(
            name=step_name, previous='Initial', timePeriod=self.time,
            massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                          BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))
        if field_outputs is None:
            field_outputs = ('S', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS',
                             'CSDMG', 'SDV', 'STATUS')
        self.model.FieldOutputRequest(
            name='F-Output-1', createStepName=step_name,
            variables=field_outputs, numIntervals=intervals)

    def create_implicit_step(self, field_outputs=None, initial_increment=0.1,
                             step_name='Loading Step'):
        self.model.StaticStep(name=step_name, previous='Initial', 
                              initialInc=initial_increment)
        if field_outputs is None:
            field_outputs = ('S', 'E', 'U', 'RF')
        self.model.FieldOutputRequest(
            name='F-Output-1', createStepName=step_name,
            variables=field_outputs)

    def create_job(self, cpus=4, subroutine_path=None, ):
        mdb.Job(
            name=self.model_name, model=self.model_name, description='',
            atTime=None, waitMinutes=0, waitHours=0, queue=None, type=ANALYSIS,
            memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
            nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',
            resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
            numDomains=cpus, activateLoadBalancing=False,
            multiprocessingMode=DEFAULT, numCpus=cpus, memory=90)

    def run_job(self):
        mdb.jobs[self.model_name].submit(consistencyChecking=OFF)
        mdb.jobs[self.model_name].waitForCompletion()

    def save_model(self):
        mdb.saveAs(pathName=self.model_name)


if __name__ == '__main__':

    # Tests
    mdl = AbaqusModel('TestModel')  # Create model
    mdl.time = 1.0  # set time attribute for testing
    mdl.create_explicit_step()
    mdl.create_job()
    mdl.save_model()
