from __future__ import division

import sys
import time
from math import sin, cos, pi, sqrt

from shelxfile import Shelxfile
from shelxfile.atoms.atoms import Atoms

"""
This module is a fork from https://github.com/des4maisons/molecule-viewer
"""


class Coordinate2D(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __mul__(self, const):
        return Coordinate2D(self.x * const, self.y * const)

    def __add__(self, coor):
        return Coordinate2D(self.x + coor.x, self.y + coor.y)

    def __sub__(self, coor):
        return self + coor * (-1)

    def __truediv__(self, const):
        return self * (1 / const)

    def __repr__(self):
        return "(%.02f, %.02f)" % (self.x, self.y)

    def __str__(self):
        return self.__repr__()


class Atom(object):
    def __init__(self, x, y, z, symbol, id):
        self.coordinate = Coordinate(x, y, z)
        self.symbol = symbol
        self.id = int(id)

    def __repr__(self):
        return str((self.symbol, self.id, self.coordinate))

    def __str__(self):
        return self.__repr__()

    def flatten(self, plane=None):
        return self.coordinate.flatten(plane)

    # Ugly! return symbol in a randomish ansi colour
    def coloured_symbol(self):
        return "\x1b[" + str(31 + (self.id % 7)) + "m" + self.symbol


class Coordinate(object):
    def __init__(self, x, y, z):
        self.x, self.y, self.z = [x, y, z]

    def __repr__(self):
        return "(%.02f, %.02f, %.02f)" % (self.x, self.y, self.z)

    def __str__(self):
        return self.__repr__()

    def __truediv__(self, const):
        return self * (1 / const)

    def __mul__(self, const):
        return Coordinate(self.x * const, self.y * const, self.z * const)

    def __add__(self, coor):
        return Coordinate(self.x + coor.x, self.y + coor.y, self.z + coor.z)

    # regular dot product
    def dot(self, vector):
        return self.x * vector.x + self.y * vector.y + self.z * vector.z

    def length(self):
        return sqrt((self.x ** 2) + (self.y ** 2) + (self.z ** 2))

    # project self onto 'onto'
    def project(self, onto):
        length = onto.length()
        scale_factor = self.dot(onto) / (length ** 2)
        return onto * scale_factor

    # cross product
    def cross(self, vector):
        a1, a2, a3 = [self.x, self.y, self.z]
        b1, b2, b3 = [vector.x, vector.y, vector.z]
        return Coordinate(a2 * b3 - a3 * b2,
                          a3 * b1 - a1 * b3,
                          a1 * b2 - a2 * b1)

    # flatten takes a 2-element list of coordinates. When interpreted as vectors,
    # these define a plane in 3-dim'l space
    def flatten(self, plane=None):
        if plane is None:  # defaults to x-y plane
            plane = [Coordinate(1, 0, 0), Coordinate(0, 1, 0)]
        # copy the list
        plane = list(plane)
        # get 2 perpendicular vectors that define the same plane
        perp = plane[0].cross(plane[1])
        plane[0] = perp.cross(plane[1])
        # make them unit vectors
        plane[0] = plane[0] / (plane[0].length())
        plane[1] = plane[1] / (plane[1].length())

        proj0 = self.project(plane[0])
        proj1 = self.project(plane[1])
        ratio0 = None
        ratio1 = None
        # the (signed) number of times the unit vector fits into the projection
        # is the 2 dimensional coordinate we want
        # so (2,0) == 2 * (1,0), 2 is what we want. don't divide by zero.
        for component in ["x", "y", "z"]:
            val = getattr(plane[0], component)
            if val != 0:
                ratio0 = getattr(proj0, component) / val
            val = getattr(plane[1], component)
            if val != 0:
                ratio1 = getattr(proj1, component) / val
        if ratio0 is None or ratio1 is None:  # if either are 0
            raise ValueError("plane defined with zero vector")
        return Coordinate2D(ratio0, ratio1)


class Molecule(object):
    def __init__(self, shx_atoms: Atoms):
        self.atoms = []
        for at in shx_atoms:
            self.atoms.append(Atom(at.x, at.y, at.z, at.name, at.resinum))

    def parse_atom(self, str):
        atype, seqnum, elt, one, x, y, z, who, cares = str.split()
        self.atoms.append(Atom(float(x), float(y), float(z), elt, int(seqnum)))

    # dimensions is 2-tuple of the number of characters on x and y axis
    def draw(self, plane=None, dimensions=None):
        if dimensions is None:
            dimensions = (80, 29)
        screen = Screen(dimensions)
        flattened = [a.flatten(plane) for a in self.atoms]
        max_extreme = Coordinate2D(max([coor.x for coor in flattened]),
                                   max([coor.y for coor in flattened]))
        min_extreme = Coordinate2D(min([coor.x for coor in flattened]),
                                   min([coor.y for coor in flattened]))
        span = max_extreme - min_extreme
        # make everything fit in our dimensions while maintaining proportions
        scale_factor = min((screen.width - 1) / span.x, (screen.height - 1) / span.y)
        extra_space = Coordinate2D(screen.width - 1, screen.height - 1) - span * scale_factor
        offset = extra_space / 2
        for atom in self.atoms:
            # shift to be in the 1st quadrant (positive coordinates) close to (0,0)
            screen_index = (atom.flatten(plane) - min_extreme) * scale_factor
            screen_index = screen_index + offset
            screen_index = (int(round(screen_index.x)), int(round(screen_index.y)))
            screen.set(screen_index, atom.coloured_symbol())
        screen.show()


class Screen(list):
    def __init__(self, widthheight):
        list.__init__(self)
        self.width = widthheight[1]
        self.height = widthheight[0]
        for i in range(self.width):
            self.append([])
            for _ in range(self.height):
                self[i].append(" ")

    def set(self, index, val):
        x = index[0]
        y = index[1]
        self[x][y] = val

    def show(self):
        for row in self:
            print("".join(row))


if __name__ == "__main__":
    shx = Shelxfile()
    shx.read_file('./tests/resources/p21c.res')
    try:
        width = sys.argv[2]
        height = sys.argv[3]
        dims = (int(width), int(height) - 1)
    except IndexError:
        dims = None
    angle = 0
    angle_incr = (2 * pi) * (1 / 40)  # 40th of a circle
    while True:
        Molecule(shx.atoms).draw(plane=[Coordinate(0, 1, 0), Coordinate(sin(angle), 0, cos(angle))], dimensions=dims)
        angle = angle + angle_incr
        time.sleep(0.5)
        print("\33[H")
