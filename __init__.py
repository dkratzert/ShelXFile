
import sys

ver = sys.version_info
if ver < (3, 4, 0):
    print("You need at least Python 3.4 to run this program!")
    sys.exit()
del sys