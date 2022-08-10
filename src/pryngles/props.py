##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#
#                                                                #
# PlanetaRY spanGLES                                             #
# The bright-side of the light-curve of (ringed) exoplanets      #
#                                                                #
##################################################################
# Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado (C) 2022  #
##################################################################
#!/usr/bin/env python
# coding: utf-8

# # Pryngles module: object properties

from pryngles import *

# ## Object properties

PROP_TYPES=["orbit","physics","optics"]
class Props():
    """
        Props is a class intended to treat the basic properties of an object 
        as an object and not as a dictionary.

        If you add another set of attributes just add the name of the attribute
        to the global variable _PROP_TYPES.
    """
    def __init__(self,**pars):
        self.__dict__.update(**pars)
    def __str__(self):
        return str(self.__dict__)


