from abaqus import *
from abaqusConstants import *
from math import pi
import mesh
import regionToolset
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\GitHub\\interlaced_model_creation\\mesoscale')
from shapely.geometry import Polygon
import sigc as sigc
import interlacedMaterials as mats
import tapePlacement as tp
from abaqus_model import AbaqusModel


class TensileModel(AbaqusModel):
    def __init__(self, model_name):
        AbaqusModel.__init__(self, model_name)

    def set_test_parameters(self, time, test_speed, output_intervals):
        '''Set test parameters'''
        self.time = time
        self.test_speed = test_speed
        self.output_intervals = output_intervals

    def set_specimen_parameters(self, t_angles, t_widths, t_spacing,
                                t_thickness, u_ratio, l_plies):
        '''
        Set interlaced specimen parameters. Note t_xxx indicates a tape
        property, u_xxx indicates an undulation property, and l_xxx indicates
        a laminate property
        '''
        self.t_angles = t_angles
        self.t_widths = [t_widths] * len(self.t_angles)
        self.t_spacing = t_spacing
        self.t_thickness = t_thickness
        self.u_ratio = u_ratio
        self.l_plies = l_plies  # number of plies in laminate
        self.sequences = (l_plies / int(len(self.t_angles))) / 2
        self.u_width = (self.t_thickness / self.u_ratio) / 2.0
        self.l_thickness = l_plies * self.t_thickness

    def create_materials(self):
        '''Create material data (uses imported matData module)'''
        E11 = 146.8
        E22 = E33 = 11.6
        nu12 = nu13 = 0.289
        nu23 = 0.298
        G12 = G13 = 6.47
        G23 = 4.38
        Xt = 2.354
        Xc = 1.102
        Yt = 0.0343
        Yc = 0.184
        Sl = 0.0827
        alpha0 = 53.0
        G1Plus = 0.1
        G1Minus = 0.1
        G2Plus = 0.00075
        G2Minus = 0.0025
        G6 = 0.0035
        density = 1.59e-06

        mats.tapeDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                        G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                        G1Minus, G2Plus, G2Minus, G6)
        mats.undulationDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12,
                              G13, G23, density, Xt, Xc, Yt, Yc, Sl, alpha0,
                              G1Plus, G1Minus, G2Plus, G2Minus, G6,
                              self.t_angles, self.t_thickness, self.u_width)
        mats.undulationDamageResin(self.model, E11, E22, E33, nu12, nu13, nu23,
                                   G12, G13, G23, density, Xt, Xc, Yt, Yc, Sl,
                                   alpha0, G1Plus, G1Minus, G2Plus, G2Minus,
                                   G6, self.t_angles, self.t_thickness,
                                   self.u_width)
        mats.resinElastic(self.model)

    def create_tape_paths(self, specimen_size):
        '''
        Uses the shapely package to define the interlaced specimen geometry.

        Args:
            specimen_size (Shapely Polygon): Polygon defining the dimensions of
                the specimen to be modeled.

        Returns:
            tape_paths (list): List of tape/tow connectivity for each tape.
            trimmed_grid (dict): Nested dictionary defining the grid of
                tape/undulation objects.
        '''
        tape_paths, part_grid = tp.laminate_creation(
            tape_angles=self.t_angles, tape_widths=self.t_widths,
            tape_spacing=self.t_spacing, specimen_size=specimen_size,
            undulation_width=self.u_width)
        return tape_paths, part_grid

    def create_3D_part(self, obj, object_id, sequence, damage=True,
                       symmetric=False):
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
            part.SectionAssignment(
                region=cell_region, sectionName=section_name,
                offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
                thicknessAssignment=FROM_SECTION)
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
                    section_name = 'Undulation-Elastic-{}'.format(section_angle)
            part.SectionAssignment(
                region=cell_region, thicknessAssignment=FROM_SECTION,
                offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
                sectionName=section_name)
        else:
            part.SectionAssignment(
                region=cell_region, sectionName='Resin-Elastic',
                offsetType=MIDDLE_SURFACE, thicknessAssignment=FROM_SECTION,
                offset=0.0, offsetField='')
        part.MaterialOrientation(
            region=cell_region, orientationType=SYSTEM, stackDirection=STACK_3,
            angle=obj.angle[-1], additionalRotationType=ROTATION_ANGLE,
            axis=AXIS_3, fieldName='', additionalRotationField='',
            localCsys=None)

        # mesh the part
        part.seedPart(size=self.mesh_size)
        part.generateMesh()

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

    def create_laminate_parts(self, part_grid, mesh_size):
        '''
        This method iterates over the part_grid dictionary and calls the
        create_3D_part method to generate parts of a laminate, one "layer" at
        a time.

        Args:
            part_grid (dict): Nested dictionary defining the grid of
                tape/undulation objects.
            mesh_size (float): Size (in mm) of the FEA mesh.

        Returns:
            None
        '''
        self.mesh_size = mesh_size
        for sequence in range(self.sequences):
            for layer_number in range(1, len(part_grid) + 1):
                for obj_id, obj in part_grid[layer_number].iteritems():
                    self.create_3D_part(obj, obj_id, sequence)

    def merge_3D_tapes(self, part_grid, tape_paths, symmetric=False):
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
                        self.merge_instances(part_grid, path_number, tape_path,
                                             sequence, mirrored=mirrored_part)
        else:
            for sequence in range(self.sequences):
                    for (path_number, tape_path) in enumerate(tape_paths):
                        self.merge_instances(part_grid, path_number, tape_path,
                                             sequence)

    def merge_instances(self, part_grid, path_number, tape_path, sequence,
                        mirrored=False):
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
                instance_name = 'Mirror-Pass3D-{}-{}'.format(sequence, path_number)
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
            pass_part.seedPart(size=self.mesh_size)
            pass_part.generateMesh()

    def get_faces_by_bounding_box(self, bounding_box):
        '''
        This method searches every instance in an assembly to determine which
        faces are contained within a bounding box.

        Args:
            bounding_box (list): A list of coodinates defining a bounding box
                in the form: (x_min, y_min, z_min, x_max, y_max, z_max)

        Returns:
            faces (list): A list of the faces within a bounding box.
        '''
        x1, y1, z1, x2, y2, z2 = bounding_box
        faces = [inst.faces.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                 for inst in self.assembly.instances.values()
                 if inst.faces.getByBoundingBox(x1, y1, z1, x2, y2, z2)]
        return faces

    def apply_constraints(self, dimensions, symmetric=False):
        '''
        This method applies constraints to the assembly.

        Args:
            dimensions (list): List of coordinates defining the exterior
                boundaries of the specimen.

        Returns:
            None
        '''

        tol = 0.001
        zMax = self.l_thickness
        zMin = 0.0
        x_min, y_min, x_max, y_max = dimensions

        # identify faces at top and bottom for coupling
        top_faces_box = (x_min - tol, y_max - tol, zMin - tol, x_max + tol,
                         y_max + tol, zMax + tol)
        top_faces = self.get_faces_by_bounding_box(top_faces_box)
        top_faces_surface = self.assembly.Surface(side1Faces=top_faces,
                                                  name='Top Faces')
        bottom_faces_box = (x_min - tol, y_min - tol, zMin - tol, x_max + tol,
                            y_min + tol, zMax + tol)
        bottom_faces = self.get_faces_by_bounding_box(bottom_faces_box)
        bottom_faces_surface = self.assembly.Surface(side1Faces=bottom_faces,
                                                     name='Bottom Faces')
        symmetric_faces_box = (x_min - tol, y_min - tol, -tol, x_max + tol,
                               y_max + tol, tol)
        symmetric_faces = self.get_faces_by_bounding_box(symmetric_faces_box)
        symmetric_faces_set = self.assembly.Set(faces=symmetric_faces,
                                                name='Symmetry Faces')

        # create reference points for coupling
        rf_point_1_id = self.assembly.ReferencePoint(
            point=(0.0, y_max + 5.0, 0.0)).id
        rf_point_1 = self.assembly.referencePoints[rf_point_1_id]
        rf_point_1_region = self.assembly.Set(referencePoints=(rf_point_1,),
                                              name='Reference Point 1')
        rf_point_2_id = self.assembly.ReferencePoint(
            point=(0.0, y_min - 5.0, 0.0)).id
        rf_point_2 = self.assembly.referencePoints[rf_point_2_id]
        rf_point_2_region = self.assembly.Set(referencePoints=(rf_point_2,),
                                              name='Reference Point 2')

        self.model.Coupling(
            name='Top Coupling', controlPoint=rf_point_1_region,
            surface=top_faces_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

        self.model.Coupling(
            name='Bottom Coupling', controlPoint=rf_point_2_region,
            surface=bottom_faces_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
            u3=ON, ur1=ON, ur2=ON, ur3=ON)

        self.model.SmoothStepAmplitude(
            name='Smoothing Amplitude', timeSpan=STEP,
            data=((0.0, 0.0), (1e-05, 1.0)))

        self.model.DisplacementBC(
            name='Bottom Surface BC', createStepName='Loading Step',
            region=rf_point_2_region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
        self.model.VelocityBC(
            name='Top Surface BC', createStepName='Loading Step',
            region=rf_point_1_region, v1=0.0, v2=self.test_speed, v3=0.0,
            vr1=0.0, vr2=0.0, vr3=0.0, amplitude='Smoothing Amplitude',
            localCsys=None, distributionType=UNIFORM, fieldName='')

        if symmetric is False:  # only create Z-symmetry if no mirrored parts
            self.model.ZsymmBC(name='Symmetry', createStepName='Loading Step',
                               region=symmetric_faces_set, localCsys=None)

    def create_interaction_properties(self):
        '''
        This method creates cohesive interaction properties.

        Args:
            None

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
        tau1 = ((9 * pi * E_33 * GIc) / (32 * Ne * self.mesh_size))**0.5
        tau2 = ((9 * pi * E_33 * GIIc) / (32 * Ne * self.mesh_size))**0.5

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

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    output_intervals = 50  # requested field output intervals
    crosshead_velocity = 0.5

    # Material parameters
    specimen_type = 'B'
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 15.0
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 8 ply laminate

    # Mesh sizes
    fineMesh = 0.125
    mediumMesh = 1.0
    coarseMesh = 3.0

    # RVE dimensions
    x_min = y_min = -(tape_widths / 2.0)
    x_max = y_max = x_min + (tape_spacing + 1) * (tape_widths)
    y_min = -75.0
    y_max = 75.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])
    dimensions = [x_min, y_min, x_max, y_max]

    # Create model
    mdl = TensileModel('Model Name')
    mdl.set_test_parameters(time, crosshead_velocity, output_intervals)
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    paths_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_laminate_parts(grid_f, mediumMesh)
    mdl.merge_3D_tapes(grid_f, paths_f)
    mdl.create_explicit_step()
    mdl.apply_constraints(dimensions)
    mdl.create_interaction_properties()
    mdl.create_interactions_explicit()
    mdl.create_job()
    mdl.save_model()
