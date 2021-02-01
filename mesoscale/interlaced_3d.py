'''
This module can be used to create 3D models of interlaced laminates.
The Interlaced3D class inherits from the InterlacedModel class. The child
class contains methods to create parts, create cohesive interactions, mesh
parts, and so forth.

(c) Rutger Kok, 25/11/2020
'''
from abaqus import *
from abaqusConstants import *
from math import pi
import mesh
import regionToolset
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\GitHub\\interlaced_model_creation\\mesoscale')
from shapely.geometry import Polygon
from interlaced_model import InterlacedModel


class Interlaced3D(InterlacedModel):
    def __init__(self, model_name, dir_name=None):
        InterlacedModel.__init__(self, model_name, dir_name)

    def create_3d_part(self, obj, object_id, sequence, damage=True,
                       symmetric=False, mesh_size=1.0):
        '''
        This function is used to create Tape parts using geometric info from
        Shapely.

        Args:
            obj (Shapely Polygon): Polygon defining the part dimensions.
            object_id (int): Grid ID of Polygon object
            sequence (): XXXXXXXX
            damage (boolean): True if the simulation includes damage i.e. is
                not purely elastic.
            symmetric (boolean): False if model uses symmetric BC, True if
                the full model is desired.

        Returns:
            None
        '''
        part_id = '{}-{}'.format(obj.layer, object_id)
        x0, y0 = obj.exterior.xy
        abaqus_coords = zip(x0, y0)
        sketch = self.model.ConstrainedSketch(
            name='Sketch {}'.format(part_id), sheetSize=200.0)
        for idx in range(len(abaqus_coords) - 1):
            sketch.Line(point1=(abaqus_coords[idx][0], abaqus_coords[idx][1]),
                        point2=(abaqus_coords[idx + 1][0],
                                abaqus_coords[idx + 1][1]))
        # extrude profile to create part
        part = self.model.Part(name='Part {}'.format(part_id),
                               dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseSolidExtrude(sketch=sketch, depth=self.t_thickness)
        cells = part.cells  # all cells within part (only 1 cell)
        cell_region = regionToolset.Region(cells=cells[:])
        if obj.object_type == 'Tape':
            if damage:
                section_name = 'Tape-Damage'
            else:
                section_name = 'Tape-Elastic'
        elif obj.object_type == 'Undulation':
            if len(obj.angle) == 1:
                if damage:
                    section_name = 'Undulation-Damage-Resin'
                else:
                    section_name = 'Undulation-Elastic-Resin'
            else:
                interface_angle = abs(obj.angle[0] - obj.angle[1])
                section_angle = (0, interface_angle)
                if damage:
                    section_name = 'Undulation-Damage-{}'.format(section_angle)
                else:
                    section_name = 'Undulation-Elastic-{}'.format(
                        section_angle)
        else:
            section_name = 'Resin-Elastic'

        part.SectionAssignment(
            region=cell_region, sectionName=section_name,
            offsetType=MIDDLE_SURFACE, thicknessAssignment=FROM_SECTION,
            offset=0.0, offsetField='')
        part.MaterialOrientation(
            region=cell_region, orientationType=SYSTEM, stackDirection=STACK_3,
            angle=obj.angle[-1], additionalRotationType=ROTATION_ANGLE,
            axis=AXIS_3, fieldName='', additionalRotationField='',
            localCsys=None)

        # mesh the part
        self.mesh_part_3d(part, mesh_size)

        # create the part instance
        instance_name = 'Instance3D-{}-{}'.format(sequence, part_id)
        sequence_offset = sequence * self.t_thickness * len(self.t_angles)
        z_offset = sequence_offset + (obj.layer - 1) * self.t_thickness
        self.assembly.Instance(name=instance_name, part=part, dependent=ON)
        self.assembly.translate(instanceList=(instance_name, ),
                                vector=(0.0, 0.0, z_offset))
        if symmetric:
            # create the symmetric part instance
            self.assembly.Instance(name='Mirror-' + instance_name, part=part,
                                   dependent=ON)
            translate_vector = (
                0.0, 0.0, self.l_thickness - z_offset - self.t_thickness)
            self.assembly.translate(
                instanceList=('Mirror-' + instance_name, ),
                vector=translate_vector)
            obj_centroid = obj.centroid.coords
            rot_point = (
                obj_centroid[0][0], obj_centroid[0][1],
                (self.l_thickness - z_offset) - self.t_thickness / 2.0)
            self.assembly.rotate(
                instanceList=('Mirror-' + instance_name, ), angle=180.0,
                axisPoint=rot_point, axisDirection=(1.0, 0.0, 0.0))

    def merge_3d_tapes(self, part_grid, tape_paths, mesh_size,
                       symmetric=False):
        '''
        This method iterates over all the tape_paths in a laminate and calls
        merge_instances to combine the separate parts into a single tape.

        Args:
            part_grid (dict): Nested dictionary defining the grid of
                tape/undulation objects.
            tape_paths (list): Nested list containing the connectivity of
                each object in each tape path.
            symmetric (boolean): False if model uses symmetric BC, True if
                the full model is desired.

        Returns:
            None
        '''
        if symmetric:
            for mirrored_part in (True, False):  # must merge mirrored parts
                for sequence in range(self.sequences):
                    for (path_number, tape_path) in enumerate(tape_paths):
                        self.merge_instances(
                            part_grid, path_number, tape_path, sequence,
                            mesh_size, mirrored=mirrored_part)
        else:
            for sequence in range(self.sequences):
                    for (path_number, tape_path) in enumerate(tape_paths):
                        self.merge_instances(
                            part_grid, path_number, tape_path, sequence,
                            mesh_size)

    def create_laminate_parts(self, part_grid, tape_paths, mesh_size,
                              symmetric=False):
        '''
        This method iterates over the part_grid dictionary and calls the
        create_3d_part method to generate parts of a laminate, one "layer" at
        a time.

        Args:
            part_grid (dict): Nested dictionary defining the grid of
                tape/undulation objects.

        Returns:
            None
        '''
        for sequence in range(self.sequences):
            for layer_number in range(1, len(part_grid) + 1):
                for obj_id, obj in part_grid[layer_number].iteritems():
                    self.create_3d_part(
                        obj, obj_id, sequence, mesh_size=mesh_size)
        self.merge_3d_tapes(part_grid, tape_paths, mesh_size, symmetric)

    def merge_instances(self, part_grid, path_number, tape_path, sequence,
                        mesh_size, mirrored=False):
        '''
        This method merges instances in the same tape path (i.e. all the
        Polygon objects that make up a single tape).

        Args:
            part_grid (dict): Nested dictionary defining the grid of
                tape/undulation objects.
            path_number (int): Number of the tape path within the current
                sequence.
            tape_path (list): List containing the connectivity of
                each object in the tape path.
            sequence (int): Number identifying the number repetitions of the
                interlacing pattern.
            mirrored (boolean): True if merging mirrored parts.

        Returns:
            None
        '''
        if mirrored:
            search_string = 'Mirror-Instance3D-{}-'.format(sequence)
        else:
            search_string = 'Instance3D-{}-'.format(sequence)

        if len(tape_path) > 1:
            layer, grid_id = tape_path[0].split('-')
            path_angle = part_grid[int(layer)][int(grid_id)].angle[-1]
            instance_list = [self.assembly.instances[search_string + str(inst)]
                             for inst in tape_path]
            if mirrored:
                instance_name = 'Mirror-Pass3D-{}-{}'.format(
                    sequence, path_number)
            else:
                instance_name = 'Pass-3D-{}-{}'.format(sequence, path_number)
            self.assembly.InstanceFromBooleanMerge(
                name=instance_name, instances=instance_list, domain=BOTH,
                keepIntersections=ON, originalInstances=DELETE, mergeNodes=ALL,
                nodeMergingTolerance=1e-06)
            pass_part = self.model.parts[instance_name]
            pass_region = regionToolset.Region(cells=pass_part.cells[:])
            pass_part.MaterialOrientation(
                region=pass_region, orientationType=SYSTEM, axis=AXIS_3,
                localCsys=None, fieldName='', stackDirection=STACK_3,
                additionalRotationType=ROTATION_ANGLE, angle=path_angle,
                additionalRotationField='')
            self.mesh_part_3d(pass_part, mesh_size)

    def mesh_part_3d(self, part, mesh_size):
        '''
        Function to mesh a part using reduced order C3D8R elements (with
        enhanced hourglass control)

        Args:
            part (Abaqus Part Instance): the part to be meshed.
            mesh_size (float): size of elements in the FE mesh.

        Returns:
            None
        '''
        elem_type_1 = mesh.ElemType(
            elemCode=C3D8R, elemLibrary=EXPLICIT, distortionControl=DEFAULT,
            kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF,
            hourglassControl=ENHANCED)
        elem_type_2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
        elem_type_3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)

        part_cells = part.Set(
            cells=part.cells[:], name='{} Cells'.format(part.name))
        part.setElementType(
            regions=part_cells,
            elemTypes=(elem_type_1, elem_type_2, elem_type_3))
        part.seedPart(size=mesh_size, deviationFactor=0.1,
                      minSizeFactor=0.1)
        part.generateMesh()

    def create_interaction_properties(self, mesh_size):
        '''
        This method creates cohesive interaction properties.

        Args:
            mesh_size (float): size of elements in FE mesh.

        Returns:
            None
        '''
        alpha = 50  # from Turon (2007)
        E_33 = 11.6
        G_12 = 6.47
        GIc = 300.0e-6
        GIIc = 800.0e-6
        Ne = 5

        # calculate cohesive zone properties according to Turon (2007)
        K1 = (alpha * E_33) / self.t_thickness
        K2 = (alpha * G_12) / self.t_thickness
        tau1 = ((9 * pi * E_33 * GIc) / (32 * Ne * mesh_size))**0.5
        tau2 = ((9 * pi * E_33 * GIIc) / (32 * Ne * mesh_size))**0.5

        self.model.ContactProperty('Tangential')
        self.model.interactionProperties['Tangential'].TangentialBehavior(
            formulation=PENALTY, directionality=ISOTROPIC,
            slipRateDependency=OFF, pressureDependency=OFF,
            temperatureDependency=OFF, dependencies=0, table=((0.15, ), ),
            shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005,
            elasticSlipStiffness=None)
        self.model.interactionProperties['Tangential'].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=ON,
            constraintEnforcementMethod=DEFAULT)
        self.model.ContactProperty('Cohesive')
        self.model.interactionProperties['Cohesive'].CohesiveBehavior(
            defaultPenalties=OFF, table=((K1, K2, K2), ))
        self.model.interactionProperties['Cohesive'].Damage(
            criterion=QUAD_TRACTION, initTable=((tau1, tau2, tau2), ),
            useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
            mixedModeType=BK, exponent=1.75, evolTable=((GIc, GIIc, GIIc), ))

    def create_interactions_explicit(self):
        '''
        This method creates cohesive interactions between all instances in an
        assembly. It uses the built-in "contactDetection" method to determine
        contacts between instances, then iterates over the identified contacts
        and adds them to a General Contact definition. This is required for
        cohesive contacts to work in Abaqus Explicit.

        Args:
            None

        Returns:
            None
        '''
        # determine contacts
        self.model.contactDetection(
            defaultType=CONTACT, interactionProperty='Cohesive',
            nameEachSurfaceFound=OFF, createUnionOfMasterSurfaces=ON,
            createUnionOfSlaveSurfaces=ON, searchDomain=MODEL,
            separationTolerance=0.0001)
        # create explicit general contact definition
        self.model.ContactExp(name='GC', createStepName='Initial')
        general_contact = self.model.interactions['GC']
        general_contact.includedPairs.setValuesInStep(
            stepName='Initial', useAllstar=ON)
        # define 'Tangential' as default behaviour
        general_contact.contactPropertyAssignments.appendInStep(
            stepName='Initial', assignments=((GLOBAL, SELF, 'Tangential'), ))
        # assign cohesive behaviour to contacting tape surfaces
        for interaction in self.model.interactions.values():
            if interaction.name == 'GC':
                continue
            else:
                master_name = interaction.master[0]
                slave_name = interaction.slave[0]
                master_surface = self.assembly.surfaces[master_name]
                slave_surface = self.assembly.surfaces[slave_name]
                general_contact.contactPropertyAssignments.appendInStep(
                    stepName='Initial',
                    assignments=((master_surface, slave_surface, 'Cohesive'), ))
                del self.model.interactions['{}'.format(interaction.name)]


if __name__ == '__main__':

    # Material parameters
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 8 ply laminate
    mesh_size = 1.0

    # RVE dimensions
    x_min = y_min = -(tape_widths / 2.0)
    x_max = y_max = x_min + (tape_spacing + 1) * (tape_widths)
    y_min = -25.0
    y_max = 25.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = Interlaced3D('TestModel')
    mdl.time = 1.0  # set time attribute for testing purposes
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    paths_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_laminate_parts(grid_f, paths_f, )
    mdl.create_interaction_properties(mesh_size)
    mdl.create_interactions_explicit()
    mdl.create_job()
    mdl.save_model()
