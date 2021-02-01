'''
This module defines the Interlaced Model class.
The class contains methods to define the parameters of the interlaced laminate,
create the tape paths and part grid from Shapely which define the geometry of
the laminate, and create the materials required.

(c) Rutger Kok, 25/11/2020
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

        mats.tapeElastic(self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                         G23, density)
        mats.tapeDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                        G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                        G1Minus, G2Plus, G2Minus, G6)
        mats.undulationDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12,
                              G13, G23, density, Xt, Xc, Yt, Yc, Sl, alpha0,
                              G1Plus, G1Minus, G2Plus, G2Minus, G6,
                              self.t_angles, self.t_thickness, self.u_width)
        mats.undulationElastic(self.model, E11, E22, E33, nu12, nu13, nu23,
                               G12, G13, G23, density, self.t_angles,
                               self.t_thickness, self.u_width)
        mats.undulationDamageResin(self.model, E11, E22, E33, nu12, nu13, nu23,
                                   G12, G13, G23, density, Xt, Xc, Yt, Yc, Sl,
                                   alpha0, G1Plus, G1Minus, G2Plus, G2Minus,
                                   G6, self.t_angles, self.t_thickness,
                                   self.u_width)
        mats.undulationElasticResin(self.model, E11, E22, E33, nu12, nu13,
                                    nu23, G12, G13, G23, density,
                                    self.t_angles, self.t_thickness,
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
        self.specimen_size = specimen_size
        tape_paths, part_grid = tp.laminate_creation(
            tape_angles=self.t_angles, tape_widths=self.t_widths,
            tape_spacing=self.t_spacing, specimen_size=specimen_size,
            undulation_width=self.u_width)
        return tape_paths, part_grid

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
    tape_widths = 15.0
    tape_spacing = 1  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 4  # symmetric 8 ply laminate

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
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    paths_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_job()
    mdl.save_model()
