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
##################################################
# PRELIMINARY INIT COMMANDS
##################################################
#--END OF TEMPLATE--#

#!/usr/bin/env python
# coding: utf-8

# # PlanetaRY spanGLES: the bright-side of the light-curve

# This is the initialization file of the `Pryngles` package.

# ## Warnings

import unittest
import warnings
import bz2
import pickle
import _pickle as cPickle
warnings.filterwarnings('ignore')

# ## Jupyter compatibilty

"""
The purpose of the get_ipython class is to provide some response in the python 
script resulting from the conversion of this notebook.

If you want to add another IPyhton function resulting from a magic command to the class, 
please verify in the resulting python script the corresponding IPython command.

For instance, the magic "%matplotlib nbagg" is converted into:

    get_ipython().magic('matplotlib nbagg',globals())

So, the routinge "magic" should be add to the get_ipython() class.        
"""
from IPython.display import HTML, Image, display
import IPython.core.autocall as autocall
from IPython import get_ipython
import sys

try:
    cfg=get_ipython().config
except AttributeError:
    def Image(url="",filename="",f=""):
        pass
    class get_ipython(object):
        def run_line_magic(self,*args):
            if "timeit" in args[0]:
                command=" ".join(args)
                replaceTimeIt(command)
        def run_cell_magic(self,x,y,z):
            pass
        def magic(self,command,scope=globals()):
            import re
            if "timeit" in command:
                replaceTimeIt(command)

#Magics can only be located from here
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

# ## PrynglesCommon
# 
# Many of the classes in Pryngles inherite methods of this common class

class PrynglesCommon(object):
    
    def save_to(self,filename,compressed=False):
        """Save object to a binary file
        
        Parameters:
            filename: string:
                Name of the file where the object will be stored.
                
            compressed: boolean, default = False:
                If True the file will be stored compressed.
        
        Notes:
            Based on https://betterprogramming.pub/load-fast-load-big-with-compressed-pickles-5f311584507e.
        """
        if compressed:
            with bz2.BZ2File(filename,"w") as f: 
                cPickle.dump(self, f)
        else:
            pikd = open(filename,"wb")
            pickle.dump(self, pikd)
            pikd.close()
            
    def load_from(self,filename,compressed=False):
        if compressed:
            pikd = bz2.BZ2File(filename,"rb")
            data = cPickle.load(pikd)
        else:
            pikd = open(filename,"rb")
            data = pickle.load(pikd)
            pikd.close()
        self.__dict__=data.__dict__
        return data
    
    def __str__(self):
        return str(self.__dict__)

# ## Miscelaneous Class

Misc_doc="""
Miscelaneous routines.

This is a set of util routines intended for a diversity of purposes.

Routines included:

    get_data(file)
""";

class Misc(object):
    def get_data(path):
        """
        Get the full path of the `datafile` which is one of the datafiles provided with the package.
        
        Parameters:
            datafile: Name of the data file, string.
            
        Return:
            Full path to package datafile in the python environment.
            
        """
        return os.path.join(ROOTDIR,'data',path);
    
    def print_df(df):
        """
        Print DataFrame.
        
        Parameters:
            df: Pandas DataFrame:
                DataFrame to print.
        """
        display(HTML(df.to_html()))
        
    def load_from(filename,compressed=False):
        
        pass
        
Misc.__doc__=Misc_doc

class Verbose(object):
    """Verbose print in the package
    
    Example:
    
        Verbose.print("Hello world") #No output
        
        Verbose.VERBOSITY=1
        Verbose.print("Hello world") #Output
    """
    VERBOSITY=0
    def print(*args):
        if Verbose.VERBOSITY:
            print(*args)

# ## Pryngles modules

from pryngles.version import *
from pryngles.consts import *
from pryngles.science import *
from pryngles.plot import *
from pryngles.props import *
from pryngles.body import *
from pryngles.sampler import *
from pryngles.spangler import *
from pryngles.star import *
from pryngles.planet import *
from pryngles.ring import *
from pryngles.observer import *
from pryngles.system import *
from pryngles.legacy import *

# ## Tests



