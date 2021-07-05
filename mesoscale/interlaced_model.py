'''
This module defines the Interlaced Model class.
The class contains methods to define the parameters of the interlaced laminate,
create the tape paths and part grid from Shapely which define the geometry of
the laminate, and create the materials required.

(c) Rutger Kok, 21/06/2021
'''
from abaqus import *
from abaqusConstants import *
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\GitHub\\interlaced_model_creation\\mesoscale')
from shapely.geometry import Polygon
import interlaced_materials as mats
import tape_placement as tp
from abaqus_model import AbaqusModel


class InterlacedModel(AbaqusModel):
    def __init__(self, model_name, dir_name=None):
        AbaqusModel.__init__(self, model_name, dir_name)

    def set_specimen_parameters(self, t_angles, t_widths, t_spacing,
                                t_thickness, u_ratio, l_plies):
        '''
        Set interlaced specimen parameters. Note t_xxx indicates a tape
        property, u_xxx indicates an undulation property, and l_xxx indicates
        a laminate property

        Args:
            t_angles (list): angles of tapes in the interlaced laminate
            t_widths (float): width of the tapes in the laminate
            t_spacing (int): number of tape width gaps between tapes placed
                in a single pass
            t_thickness (float): thickness of the tapes in the laminate
            u_ratio (float): ratio of the height of an undulation to its length
            l_plies (int): number of plies in the laminate.

        Returns:
            None
        '''
        self.t_angles = t_angles
        self.t_widths = [t_widths] * len(self.t_angles)
        self.t_spacing = t_spacing
        self.t_thickness = t_thickness
        self.u_ratio = u_ratio
        self.l_plies = l_plies  # number of plies in laminate
        self.sequences = int((l_plies / int(len(self.t_angles))) / 2)
        self.u_width = (self.t_thickness / self.u_ratio) / 2.0
        self.l_thickness = l_plies * self.t_thickness

    def create_materials(self):
        '''
        Create materials required for interlaced models. Makes use of the
        interlaced_materials module.

        Args:
            None

        Returns:
            None
        '''
        # Material input data
        E11 = 116.6
        E22 = E33 = 7.231
        nu12 = nu13 = 0.339
        nu23 = 0.374
        G12 = G13 = 3.268
        G23 = 2.632
        Xt = 2.180
        Xpo = 0.218  # estimate based on Maimi
        Xc = 0.811
        Yt = 0.131
        Yc = 0.185
        Sl = 0.122
        St = 0.070
        etaL = 0.5
        alpha0 = 53.0
        GL1Plus = 0.1034  # Tan
        GE1Plus = 0.0296  # Estimate from Tan and Maimi (same ratio)
        G1Minus = 0.095  # Furtado
        G2Plus = 0.00038
        G6 = 0.00162
        density = 1.59e-06

        mats.tape_elastic(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            density)
        mats.tape_damage(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23, Xt,
            Xpo, Xc, Yt, Yc, Sl, St, alpha0, etaL, GL1Plus, GE1Plus, G1Minus,
            G2Plus, G6, density)
        mats.undulation_damage(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23, Xt,
            Xpo, Xc, Yt, Yc, Sl, St, alpha0, etaL, GL1Plus, GE1Plus, G1Minus,
            G2Plus, G6, density, self.t_angles, self.t_thickness, self.u_width)
        mats.undulation_elastic(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            density, self.t_angles, self.t_thickness, self.u_width)
        mats.resin_elastic(
            self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23,
            density)
        mats.create_interaction_properties(
            self.model, E33, G12, G2Plus, G6, self.mesh_size, self.t_thickness)

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
        self.specimen_size = specimen_size
        x, y = zip(*self.specimen_size.exterior.coords)
        self.x_max = max(x)  # determine max x-value
        self.y_max = max(y)
        self.x_min = min(x)
        self.y_min = min(y)
        tape_paths, tape_angles, part_grid = tp.laminate_creation(
            tape_angles=self.t_angles, tape_widths=self.t_widths,
            tape_spacing=self.t_spacing, specimen_size=specimen_size,
            undulation_width=self.u_width)
        return tape_paths, tape_angles, part_grid

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

    def get_edges_by_bounding_box(self, bounding_box):
        '''
        This method searches every instance in an assembly to determine which
        edges are contained within a bounding box.

        Args:
            bounding_box (list): A list of coodinates defining a bounding box
                in the form: (x_min, y_min, z_min, x_max, y_max, z_max)

        Returns:
            edges (list): A list of the edges within a bounding box.
        '''
        x1, y1, z1, x2, y2, z2 = bounding_box
        edges = [inst.edges.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                 for inst in self.assembly.instances.values()
                 if inst.edges.getByBoundingBox(x1, y1, z1, x2, y2, z2)]
        return edges

    def get_cells_by_bounding_box(self, bounding_box):
        '''
        This method searches every instance in an assembly to determine which
        cells are contained within a bounding box.

        Args:
            bounding_box (list): A list of coodinates defining a bounding box
                in the form: (x_min, y_min, z_min, x_max, y_max, z_max)

        Returns:
            cells (list): A list of the cells within a bounding box.
        '''
        x1, y1, z1, x2, y2, z2 = bounding_box
        cells = [inst.cells.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                 for inst in self.assembly.instances.values()
                 if inst.cells.getByBoundingBox(x1, y1, z1, x2, y2, z2)]
        return cells


if __name__ == '__main__':
    # Material parameters
    tape_angles = (0, 90)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 4 ply laminate

    # RVE dimensions
    x_min = y_min = -(tape_widths / 2.0)
    x_max = y_max = x_min + (tape_spacing + 1) * (tape_widths)
    y_min = -75.0
    y_max = 75.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = InterlacedModel('TestModel')
    mdl.time = 1.0  # set time attribute for testing purposes
    mdl.mesh_size = 0.5  # set mesh size for testing purposes
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    paths_f, angles_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_job()
    mdl.save_model()
