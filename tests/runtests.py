import doctest

import dsrmath
import misc
from refine import shx_refine
from shelxfile import shelx, cards, atoms

modules = [shelx, cards, misc, atoms, shx_refine] 


def run_tests():
    for name in modules:
        print("Testing {} ...".format(name.__name__))
        failed, attempted = doctest.testmod(name)  # , verbose=True)
        if failed == 0:
            print('passed all {} tests in {}!'.format(attempted, name.__name__))
        else:
            msg = '!!!!!!!!!!!!!!!! {} of {} tests failed in {}  !!!!!!!!!!!!!!!!'.format(failed,
                                                                                          attempted,
                                                                                          name.__name__)

if __name__ == '__main__':
    run_tests()