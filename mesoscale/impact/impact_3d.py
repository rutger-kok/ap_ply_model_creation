'''
Module used to create a model of a Drop Weight Tower test on an interlaced
laminate. Uses the Interlaced3D class to create the interlaced laminate
geometry.

(c) Rutger Kok, 27/11/2020
'''
from abaqus import *
from abaqusConstants import *
import mesh
from sys import path
githubPath = 'C:\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\mesoscale')
path.append(githubPath + '\\interlaced_model_creation\\mesoscale\\impact')
from shapely.geometry import Polygon
from impact_model import ImpactModel


class DWTModel(ImpactModel):
    def __init__(self, model_name, dir_name=None):
        ImpactModel.__init__(self, model_name, dir_name)

    def create_shell_part(self, x_interior, y_interior, x_exterior, y_exterior,
                          part_grid):
        '''
        This function is used to create a (single!) 2D shell Tape part using
        geometric info from Shapely.
        '''
        # Create sketch
        sketch = self.model.ConstrainedSketch(
            name='2D Part Sketch', sheetSize=200.0)
        sketch.rectangle(point1=(-x_exterior, -y_exterior),
                         point2=(x_exterior, y_exterior))
        sketch.rectangle(point1=(-x_interior, -y_interior),
                         point2=(x_interior, y_interior))
        part = self.model.Part(
            name='Shell Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseShell(sketch=sketch)

        self.partition_part(part)
        self.assign_properties(part, part_grid)
        self.mesh_part_2d(part, self.medium_mesh)

        instance_name = 'Shell Part Instance'
        self.assembly.Instance(name=instance_name, part=part, dependent=ON)
        self.assembly.rotate(
            instanceList=(instance_name, ), angle=-90.0,
            axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0))
        self.shell_edge_point = (x_exterior, 0.0, 0.0)

    def rotate_3d_instances(self):
        non_3d_instances = ['Impactor Instance', 'Bottom Plate Instance',
                            'Clamp Instance 1', 'Clamp Instance 2',
                            'Clamp Instance 3', 'Clamp Instance 4',
                            'Shell Part Instance']
        instances_3d = [inst.name for inst in self.assembly.instances.values()
                        if inst.name not in non_3d_instances]
        self.assembly.rotate(instanceList=instances_3d, angle=-90.0,
                             axisPoint=(0.0, 0.0, 0.0),
                             axisDirection=(1.0, 0.0, 0.0))

    def tie_shell_to_3d(self, x, y):
        x_min = x  # in rotated frame!
        x_max = -x
        z_min = y
        z_max = -y
        y_min = 0.0
        y_max = self.l_thickness
        tol = 0.001

        shell_instance = self.assembly.instances['Shell Part Instance']
        shell_edge = shell_instance.edges.findAt(
            ((x_min, 0.0, 0.0), ), ((x_max, 0.0, 0.0), ),
            ((0.0, 0.0, z_min), ), ((0.0, 0.0, z_max), ))
        tie_region_shell = self.assembly.Surface(
            side1Edges=shell_edge, name='Tie Edges Shell')

        side_1_bbox = [x_max - tol, y_min - tol, z_min - tol, x_max + tol,
                       y_max + tol, z_max + tol]
        side_1_faces = self.get_faces_by_bounding_box(side_1_bbox)
        side_2_bbox = [x_min - tol, y_min - tol, z_max - tol, x_max + tol,
                       y_max + tol, z_max + tol]
        side_2_faces = self.get_faces_by_bounding_box(side_2_bbox)
        side_3_bbox = [x_min - tol, y_min - tol, z_min - tol, x_min + tol,
                       y_max + tol, z_max + tol]
        side_3_faces = self.get_faces_by_bounding_box(side_3_bbox)
        side_4_bbox = [x_min - tol, y_min - tol, z_min - tol, x_max + tol,
                       y_max + tol, z_min + tol]
        side_4_faces = self.get_faces_by_bounding_box(side_4_bbox)
        all_sides = side_1_faces + side_2_faces + side_3_faces + side_4_faces
        tie_region_3d = self.assembly.Surface(
            side1Faces=all_sides, name='Tie Faces 3D')
        self.model.Tie(
            name='Shell-3D-Tie', master=tie_region_3d, slave=tie_region_shell,
            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON,
            thickness=ON)

    def create_interactions_explicit(self):
        # determine contacts
        non_cohesive_instances = ['Impactor Instance', 'Bottom Plate Instance',
                                  'Clamp Instance 1', 'Clamp Instance 2',
                                  'Clamp Instance 3', 'Clamp Instance 4',
                                  'Shell Part Instance']
        cohesive_instances = []
        for inst in self.assembly.instances.values():
            if inst.name not in non_cohesive_instances:
                cohesive_instances.append(inst.name)
        self.model.contactDetection(
            defaultType=CONTACT, interactionProperty='Cohesive',
            nameEachSurfaceFound=OFF, createUnionOfMasterSurfaces=ON,
            createUnionOfSlaveSurfaces=ON, searchDomain=cohesive_instances,
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
    energy = 30.0  # impact energy to simulate

    # Material parameters
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 8 ply laminate

    # Mesh sizes
    fine_mesh = 0.25
    medium_mesh = 1.0
    coarse_mesh = 3.0

    # Define specimen dimensions
    # Fine mesh region
    xMin_f = -25.0
    xMax_f = -xMin_f
    yMin_f = -25.0
    yMax_f = -yMin_f
    fine_region = Polygon([(xMin_f, yMin_f), (xMax_f, yMin_f),
                           (xMax_f, yMax_f), (xMin_f, yMax_f)])

    # Shell region
    x_interior = -25.0
    x_exterior = -50.0
    y_interior = -25.0
    y_exterior = 75.0
    exterior = [(x_exterior, y_exterior), (x_exterior, -y_exterior),
                (-x_exterior, -y_exterior), (-x_exterior, y_exterior),
                (x_exterior, y_exterior)]
    interior = [(x_interior, y_interior), (x_interior, -y_interior),
                (-x_interior, -y_interior), (-x_interior, y_interior),
                (x_interior, y_interior)][::-1]
    shell_region = Polygon(exterior, [interior])

    # Create model
    mdl = DWTModel('DWT_Test')
    mdl.set_test_parameters(time, output_intervals, energy, coarse_mesh,
                            medium_mesh, fine_mesh)
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    mdl.create_impactor_part()
    mdl.create_clamp_part()
    mdl.create_bottom_plate()

    paths_sh, grid_sh = mdl.create_tape_paths(shell_region)
    mdl.create_shell_part(x_interior, y_interior, x_exterior, y_exterior,
                          grid_sh)
    paths_f, grid_f = mdl.create_tape_paths(fine_region)
    mdl.create_laminate_parts(grid_f, paths_f, coarse_mesh)
    mdl.create_explicit_step()
    mdl.create_assembly()
    mdl.rotate_3d_instances()
    mdl.constrain_b_plate()
    mdl.constrain_clamps()
    mdl.constrain_impactor()
    mdl.tie_shell_to_3d(x_interior, y_interior)
    mdl.create_interaction_properties(coarse_mesh)
    # mdl.create_interactions_explicit()
    mdl.create_job()
    mdl.save_model()
