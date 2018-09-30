
from .shelxfile.shelx import __version__
from .shelxfile.shelx import __doc__
from .shelxfile.shelx import *
from .shelxfile.atoms import *
from .shelxfile.misc import *
from .shelxfile.cards import *
from .shelxfile.dsrmath import *
from .refine.shx_refine import *

__all__ = [
    "ShelXFile",
    "Atoms",
    "Atom",
    "frac_to_cart",
    "cart_to_frac",
    "SDM",
    "elements"]

if __name__ == "__main__":
    ver = sys.version_info
    if ver < (3, 4, 0):
        print("You need at least Python 3.4 to run this program!")
        sys.exit()
    del sys


