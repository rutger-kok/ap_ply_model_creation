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
sys.path.append('C:\\Github\\interlaced_model_creation')
import numpy as np
from abaqus import *
from abaqusConstants import *
from shapely.geometry import Polygon
from ap_ply_3d import AP_PLY_3D


class AP_PLY_Tensile_3D(AP_PLY_3D):
    def __init__(self, model_name, dir_name=None):
        AP_PLY_3D.__init__(self, model_name, dir_name)

    def set_test_parameters(self, time, test_speed, output_intervals,
                            solid_mesh, symmetry_mode):
        '''Set test parameters'''
        self.time = time  # simulation time in ms
        self.test_speed = test_speed  # crosshead velocity in mm/s
        self.output_intervals = output_intervals  # number of output intervals
        self.solid_mesh = solid_mesh  # mesh size in mm
        self.symmetry_mode = symmetry_mode  # 'BC' or 'GEOM'

    def apply_constraints(self):
        '''
        This method applies constraints to the assembly.

        Args:
            None

        Returns:
            None
        '''
        tol = 0.001
        # identify faces at top and bottom for coupling
        right_faces_bbox = (
            50.0 - tol, -32.5 - tol, -self.l_thickness - tol,
            50.0 + tol, 32.5 + tol, self.l_thickness + tol)
        right_faces = self.get_faces_by_bounding_box(right_faces_bbox)
        right_faces_surface = self.assembly.Surface(
            side1Faces=right_faces, name='Right Faces')
        left_faces_bbox = (
            -50.0 - tol, -32.5 - tol, -self.l_thickness - tol,
            -50.0 + tol, 32.5 + tol, self.l_thickness + tol)
        left_faces = self.get_faces_by_bounding_box(left_faces_bbox)
        left_faces_surface = self.assembly.Surface(
            side1Faces=left_faces, name='Left Faces')

        # create reference points for coupling
        rf_point_1_id = self.assembly.ReferencePoint(
            point=(60.0, 0.0, 0.410)).id
        rf_point_1 = self.assembly.referencePoints[rf_point_1_id]
        rf_point_1_region = self.assembly.Set(
            referencePoints=(rf_point_1,), name='Reference Point 1')
        rf_point_2_id = self.assembly.ReferencePoint(
            point=(-60.0, 0.0, 0.410)).id
        rf_point_2 = self.assembly.referencePoints[rf_point_2_id]
        rf_point_2_region = self.assembly.Set(
            referencePoints=(rf_point_2,), name='Reference Point 2')

        self.model.Coupling(
            name='Right Coupling', controlPoint=rf_point_1_region,
            surface=right_faces_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

        self.model.Coupling(
            name='Left Coupling', controlPoint=rf_point_2_region,
            surface=left_faces_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
            u3=ON, ur1=ON, ur2=ON, ur3=ON)

        # resultant velocity BC
        # self.model.DisplacementBC(
        #     name='Left Surface BC', createStepName='Loading Step',
        #     region=rf_point_2_region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
        #     ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
        #     fieldName='', localCsys=None)
        # self.model.TabularAmplitude(
        #     name='V_in', timeSpan=STEP, smooth=SOLVER_DEFAULT,
        #     data=(
        #         (0.0, 0.0), (0.05, 0.47438150020065),
        #         (0.1, 1.46392861552735), (0.15, 2.68617243159958),
        #         (0.2, 3.41520540580244), (0.25, 3.82353305499804),
        #         (0.3, 3.99067287265618), (0.35, 4.31825736080335),
        #         (0.4, 4.47221614475103), (0.45, 4.3580065689109),
        #         (0.5, 4.39284784454535), (0.55, 4.11785927154367),
        #         (0.6, 4.10578819146008), (0.65, 3.99013246214402)))

        # velocity BC both sides
        self.model.TabularAmplitude(
            name='V_in', timeSpan=STEP, smooth=SOLVER_DEFAULT,
            data=[
                (0.0, 0.00039366), (0.05, 0.47508999999999996),
                (0.1, 1.4763), (0.15, 2.736), (0.2, 3.5184),
                (0.25, 4.0149), (0.3, 4.2783), (0.35, 4.6886),
                (0.4, 4.9355), (0.45, 4.9275), (0.5, 5.038),
                (0.55, 4.8216), (0.6, 4.852), (0.65, 4.7943)])
        self.model.TabularAmplitude(
            name='V_out', timeSpan=STEP, smooth=SOLVER_DEFAULT,
            data=[
                (0.0, 0.003865), (0.05, 0.00070438), (0.1, 0.012326),
                (0.15, 0.049843), (0.2, 0.10318), (0.25, 0.19138),
                (0.3, 0.28764), (0.35, 0.37035), (0.4, 0.46328),
                (0.45, 0.56953), (0.5, 0.64518), (0.55, 0.70373),
                (0.6, 0.74622), (0.65, 0.8042)])

        self.model.VelocityBC(
            name='Right Surface BC', createStepName='Loading Step',
            region=rf_point_1_region, v1=1.0, v2=0.0, v3=0.0,
            vr1=0.0, vr2=0.0, vr3=0.0, amplitude='V_in',
            localCsys=None, distributionType=UNIFORM, fieldName='')

        self.model.VelocityBC(
            name='Left Surface BC', createStepName='Loading Step',
            region=rf_point_2_region, v1=1.0, v2=0.0, v3=0.0,
            vr1=0.0, vr2=0.0, vr3=0.0, amplitude='V_out',
            localCsys=None, distributionType=UNIFORM, fieldName='')

        if self.symmetry_mode == 'BC':
            symmetric_faces_bbox = (
                -60.0, -40.0, -tol, 60.0, 40.0, tol)
            symmetric_faces = self.get_faces_by_bounding_box(
                symmetric_faces_bbox)
            symmetric_faces_set = self.assembly.Set(
                faces=symmetric_faces, name='Symmetry Faces')
            self.model.ZsymmBC(
                name='Symmetry', createStepName='Loading Step',
                region=symmetric_faces_set, localCsys=None)

    def specimen_extraction(self, Lo, wo, D, wg, Lg):
        x, y = zip(*self.specimen_size.exterior.coords)
        self.x_max = max(x)
        self.x_min = min(x)
        self.y_max = max(y)
        self.y_min = min(y)
        x0 = ((self.x_max + self.x_min) / 2.0) + self.x_shift
        y0 = ((self.y_max + self.y_min) / 2.0) + self.y_shift
        partSketch = self.model.ConstrainedSketch(
            name='Profile', sheetSize=200.0)
        g = partSketch.geometry
        # leftmost edge
        l1 = partSketch.Line(
            point1=(x0 - Lo / 2.0, y0),
            point2=(x0 - Lo / 2.0, y0 + wo / 2.0))
        # top left
        l2 = partSketch.Line(
            point1=(x0 - Lo / 2.0, y0 + wo / 2.0),
            point2=(x0 - D / 2.0, y0 + wo / 2.0))
        # top right
        l3 = partSketch.Line(
            point1=(x0 + D / 2.0, y0 + wo / 2.0),
            point2=(x0 + Lo / 2.0, y0 + wo / 2.0))
        # rightmost edge
        l4 = partSketch.Line(
            point1=(x0 + Lo / 2.0, y0),
            point2=(x0 + Lo / 2.0, y0 + wo / 2.0))
        # gauge length line
        l5 = partSketch.Line(
            point1=(x0 - Lg / 2.0, y0 + wg / 2.0),
            point2=(x0 + Lg / 2.0, y0 + wg / 2.0))
        # shoulder
        arc1 = partSketch.ArcByStartEndTangent(
            point1=(x0 - Lg / 2.0, y0 + wg / 2.0),
            point2=(x0 - D / 2.0, y0 + wo / 2.0),
            entity=l5)
        arc2 = partSketch.ArcByStartEndTangent(
            point1=(x0 + Lg / 2.0, y0 + wg / 2.0),
            point2=(x0 + D / 2.0, y0 + wo / 2.0),
            entity=l5)

        partSketch.ConstructionLine(
            point1=(0.0, y0), point2=(100.0, y0))
        partSketch.copyMirror(
            mirrorLine=g.findAt((50.0, y0)),
            objectList=(g.findAt(arc1.pointOn),
                        g.findAt(arc2.pointOn),
                        g.findAt(l1.pointOn),
                        g.findAt(l2.pointOn),
                        g.findAt(l3.pointOn),
                        g.findAt(l4.pointOn),
                        g.findAt(l5.pointOn), ))

        for part in self.model.parts.values():
            vertex_array = np.asarray([v.pointOn[0] for v in part.vertices])
            xmax, ymax, zmax = vertex_array.max(axis=0)
            xmin, ymin, zmin = vertex_array.min(axis=0)
            sketch_plane = part.faces.findAt(
                (xmax - 0.0001, ymin + 0.0001, zmax))
            sketch_edge = part.edges.findAt((xmax, ymin + 0.0001, zmax))
            t = part.MakeSketchTransform(
                sketchPlane=sketch_plane,
                sketchUpEdge=sketch_edge,
                sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
                origin=(x0, y0, zmax))
            s = self.model.ConstrainedSketch(
                name='Sweep', sheetSize=200.0, gridSpacing=20.0, transform=t)
            s.retrieveSketch(sketch=self.model.sketches['Profile'])
            s.rectangle(
                point1=(xmin - 150.0, ymin - 150.0),
                point2=(xmax + 150.0, ymax + 150.0))
            part.CutExtrude(
                sketchPlane=sketch_plane,
                sketchUpEdge=sketch_edge,
                sketchPlaneSide=SIDE1,
                sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=OFF)


if __name__ == '__main__':
    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    output_intervals = 100  # requested field output intervals
    test_speed = 1.0  # [mm/ms]
    symmetry_mode = 'BC'
    solid_mesh = 1.5

    # Material parameters
    # tape_angles = (0, 45, 90, -45)   # define angles of tapes
    tape_angles = (0, 90)   # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.205  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.0683333  # ratio of undulation amplitude to length
    number_of_plies = 8  # symmetric 4 ply laminate

    # RVE dimensions
    x_min = -50.0
    x_max = 50.0
    y_min = -32.5
    y_max = 32.5
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])
    # Create model
    mdl = AP_PLY_Tensile_3D('XP_ELSA_vinvout')
    mdl.set_test_parameters(time, test_speed, output_intervals,
                            solid_mesh, symmetry_mode)
    mdl.create_explicit_step()
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies, x_shift=0.0, y_shift=0.0)
    mdl.create_materials()
    paths_f, angles_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_parts(grid_f, mesh=False)
    mdl.create_sequences()
    mdl.specimen_extraction(250.0, 65.0, 100.0, 40.0, 40.0)
    mdl.create_cohesive_interactions(tie=True)
    mdl.apply_constraints()
    mdl.create_job()
    mdl.save_model()
