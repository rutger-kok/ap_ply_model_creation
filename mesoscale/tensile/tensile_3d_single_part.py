from abaqus import *
from abaqusConstants import *
from math import pi
import mesh
import regionToolset
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\GitHub\\interlaced_model_creation\\mesoscale')
from shapely.geometry import Polygon
from interlaced_2d import Interlaced2D
from interlaced_3d import Interlaced3D


class TensileModel(Interlaced2D, Interlaced3D):
    '''
    Inherits partition_part from Interlaced2D and mesh_part_3d from
    Interlaced3D.
    '''
    def __init__(self, model_name, dir_name=None):
        Interlaced2D.__init__(self, model_name, dir_name)

    def set_test_parameters(self, time, test_speed, output_intervals,
                            mesh_size=0.5):
        '''Set test parameters'''
        self.time = time
        self.test_speed = test_speed
        self.output_intervals = output_intervals
        self.mesh_size = mesh_size

    def create_part(self):

        x, y = zip(*self.specimen_size.exterior.coords)
        x_max = max(x)
        y_max = max(y)
        x_min = min(x)
        y_min = min(y)
        z_max = self.l_thickness

        sketch = self.model.ConstrainedSketch(
            name='Specimen Sketch', sheetSize=200.0)
        sketch.rectangle(point1=(x_min, y_min), point2=(x_max, y_max))
        self.specimen_part = self.model.Part(
            name='Specimen', dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        self.specimen_part.BaseSolidExtrude(sketch=sketch, depth=z_max)
        specimen_instance = self.assembly.Instance(
            name='Specimen', part=self.specimen_part, dependent=ON)

    def partition_through_thickness(self):
        datum_plane_ids = []
        for i in range(1, self.l_plies):
            datum_plane_ids.append(
                self.specimen_part.DatumPlaneByPrincipalPlane(
                    principalPlane=XYPLANE, offset=self.t_thickness * i).id)
            specimen_cells = self.specimen_part.cells
            self.specimen_part.PartitionCellByDatumPlane(
                datumPlane=self.specimen_part.datums[datum_plane_ids[i - 1]],
                cells=specimen_cells)

    def assign_properties(self, part_grid, damage=True):
        for layer, layer_grid in part_grid.iteritems():
            for obj_id, obj in layer_grid.iteritems():
                obj_centroid = obj.centroid.coords
                z_coord = (2 * layer - 1) / 2.0 * self.t_thickness
                z_coords = [z_coord + n * (self.t_thickness * len(self.t_angles))
                            for n in range(0, self.sequences + 1)]
                for i, z in enumerate(z_coords):
                    if i == 0:
                        selected_cells = self.specimen_part.cells.findAt(
                            ((obj_centroid[0][0], obj_centroid[0][1], z), ))
                    else:
                        selected_cells += self.specimen_part.cells.findAt(
                            ((obj_centroid[0][0], obj_centroid[0][1], z), ))
                obj_region = regionToolset.Region(cells=selected_cells)
                obj_angle = obj.angle[-1]
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
                            section_name = 'Undulation-Damage-{}'.format(
                                section_angle)
                        else:
                            section_name = 'Undulation-Elastic-{}'.format(
                                section_angle)
                else:
                    section_name = 'Resin-Elastic'
                self.specimen_part.SectionAssignment(
                    region=obj_region, sectionName=section_name,
                    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
                    thicknessAssignment=FROM_SECTION)
                self.specimen_part.MaterialOrientation(
                    region=obj_region, orientationType=SYSTEM, axis=AXIS_3,
                    localCsys=None, fieldName='', stackDirection=STACK_3,
                    additionalRotationType=ROTATION_ANGLE,
                    additionalRotationField='', angle=obj_angle, )
        self.mesh_part_3d(self.specimen_part, self.mesh_size)

    def apply_constraints(self, symmetric=False):
        '''
        This method applies constraints to the assembly.

        Args:
            symmetric (boolean): False if the assembly does not contain
                mirrored instances (i.e. symmetry is enforced through a
                constraint rather than through geometric symmetry)

        Returns:
            None
        '''
        tol = 0.001
        x, y = zip(*self.specimen_size.exterior.coords)
        x_max = max(x)  # determine max x-value for offsets
        y_max = max(y)
        x_min = min(x)  # determine max x-value for offsets
        y_min = min(y)
        z_min = 0.0
        z_max = self.l_thickness

        # identify faces at top and bottom for coupling
        top_faces_box = (x_min - tol, y_max - tol, z_min - tol, x_max + tol,
                         y_max + tol, z_max + tol)
        top_faces = self.get_faces_by_bounding_box(top_faces_box)
        top_faces_surface = self.assembly.Surface(side1Faces=top_faces,
                                                  name='Top Faces')
        bottom_faces_box = (x_min - tol, y_min - tol, z_min - tol, x_max + tol,
                            y_min + tol, z_max + tol)
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


if __name__ == '__main__':

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    output_intervals = 50  # requested field output intervals
    crosshead_velocity = 0.5

    # Material parameters
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 8 ply laminate

    # Mesh sizes
    mesh_size = 2.0

    # RVE dimensions
    x_min = -(tape_widths / 2.0)
    x_max = x_min + (tape_spacing + 1) * (tape_widths)
    y_min = -25.0
    y_max = 25.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = TensileModel('TestModel')
    mdl.set_test_parameters(time, crosshead_velocity, output_intervals)
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    paths_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_part()
    mdl.partition_through_thickness()
    mdl.partition_part(face=False)
    mdl.create_explicit_step()
    mdl.apply_constraints()
    mdl.create_job()
    mdl.save_model()
