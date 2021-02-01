
from abaqusConstants import *
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\GitHub\\interlaced_model_creation\\mesoscale\\tensile')
from shapely.geometry import Polygon
from tensile_3d_single_part import TensileModel


class ELSA(TensileModel):
    def __init__(self, model_name, dir_name=None):
        TensileModel.__init__(self, model_name, dir_name)

    def trim_specimen(self):
        sketch_plane = self.specimen_part.faces.findAt(
            (0.0, 0.0, self.l_thickness), )
        sketch_edge = self.specimen_part.edges.findAt(
            (self.t_widths[0] * 0.5 - self.u_width, 0.0, self.l_thickness), )
        sketch_transform = self.specimen_part.MakeSketchTransform(
            sketchPlane=sketch_plane, sketchUpEdge=sketch_edge,
            sketchPlaneSide=SIDE1, origin=(-32.5, -125.0, self.l_thickness))
        sketch = self.model.ConstrainedSketch(
            name='Specimen Trim', sheetSize=20.0, gridSpacing=0.5,
            transform=sketch_transform)
        lines = [((0, 75), (0, 175)), ((65, 75), (65, 175))]
        for line in lines:
            sketch.Line(point1=line[0], point2=line[1])
        gauge_line_1 = sketch.Line(point1=(12.5, 105), point2=(12.5, 145))
        gauge_line_2 = sketch.Line(point1=(52.5, 105), point2=(52.5, 145))
        arc_1 = sketch.ArcByStartEndTangent(
            point1=(12.5, 105), point2=(0, 75), entity=gauge_line_1)
        arc_2 = sketch.ArcByStartEndTangent(
            point1=(12.5, 145), point2=(0, 175), entity=gauge_line_1)
        arc_3 = sketch.ArcByStartEndTangent(
            point1=(52.5, 105), point2=(65, 75), entity=gauge_line_2)
        arc_4 = sketch.ArcByStartEndTangent(
            point1=(52.5, 145), point2=(65, 175), entity=gauge_line_2)
        self.specimen_part.CutExtrude(sketchPlane=sketch_plane,
        sketchUpEdge=sketch_edge, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=sketch, flipExtrudeDirection=OFF)

        self.mesh_part_3d(self.specimen_part, self.mesh_size)  # remesh

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
        # origin is at center of gauge section
        datum_plane_1_id = self.specimen_part.DatumPlaneByPrincipalPlane(
            principalPlane=XZPLANE, offset=50.0).id
        datum_plane_2_id = self.specimen_part.DatumPlaneByPrincipalPlane(
            principalPlane=XZPLANE, offset=-50.0).id
        self.specimen_part.PartitionCellByDatumPlane(
            datumPlane=self.specimen_part.datums[datum_plane_1_id],
            cells=self.specimen_part.cells)
        self.specimen_part.PartitionCellByDatumPlane(
            datumPlane=self.specimen_part.datums[datum_plane_2_id],
            cells=self.specimen_part.cells)

        tol = 0.001
        z_min = 0.0
        z_max = self.l_thickness

        # identify faces at top and bottom for coupling
        top_cells_box = (-32.5 - tol, 50.0 - tol, z_min - tol, 32.5 + tol,
                         125.0 + tol, z_max + tol)
        top_cells = self.get_cells_by_bounding_box(top_cells_box)
        top_cells_surface = self.assembly.Set(cells=top_cells,
                                              name='Top Cells')
        bottom_cells_box = (-32.5 - tol, -125.0 - tol, z_min - tol, 32.5 + tol,
                            -50.0 + tol, z_max + tol)
        bottom_cells = self.get_cells_by_bounding_box(bottom_cells_box)
        bottom_cells_surface = self.assembly.Set(cells=bottom_cells,
                                                 name='Bottom Cells')
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
            surface=top_cells_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

        self.model.Coupling(
            name='Bottom Coupling', controlPoint=rf_point_2_region,
            surface=bottom_cells_surface, influenceRadius=WHOLE_SURFACE,
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
    tape_widths = 25.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 8 ply laminate

    # Mesh sizes
    mesh_size = 5.0

    # RVE dimensions
    x_min = -32.5
    x_max = 32.5
    y_min = -125.0
    y_max = 125.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = ELSA('TestModel')
    mdl.set_test_parameters(time, crosshead_velocity, output_intervals,
                            mesh_size=mesh_size)
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    paths_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_part()
    mdl.partition_through_thickness()
    mdl.partition_part(face=False)
    mdl.assign_properties(grid_f)
    mdl.trim_specimen()
    mdl.create_explicit_step()
    mdl.apply_constraints()
    mdl.create_job()
    mdl.save_model()
