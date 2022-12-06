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
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# External required packages
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import unittest
import warnings
import dill
import inspect
import sigfig
from copy import deepcopy
import sys
from collections import OrderedDict as odict
warnings.filterwarnings('ignore')

#JupDev: Jupyter compatibility
from IPython.display import HTML, Image, display
import IPython.core.autocall as autocall
from IPython import get_ipython

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Stand alone code of the module
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

"""
The purpose of the get_ipython class is to provide some response in the python 
script resulting from the conversion of this notebook.

If you want to add another IPyhton function resulting from a magic command to the class, 
please verify in the resulting python script the corresponding IPython command.

For instance, the magic "%matplotlib nbagg" is converted into:

    get_ipython().magic('matplotlib nbagg',globals())

So, the method "magic" should be add to the get_ipython() class.        
"""
try:
    cfg=get_ipython().config
except AttributeError:
    def Image(url="",filename="",f=""):
        pass
    class get_ipython(object):
        def run_line_magic(self,*args):
            pass
        def run_cell_magic(self,x,y,z):
            pass
        def magic(self,command,scope=globals()):
            pass

#Magics can only be located starting from here
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

#Verbosity levels: see help(Verbose)
VERB_NONE=0
VERB_SIMPLE=1
VERB_SYSTEM=2
VERB_VERIFY=3
VERB_DEEP=4
VERB_ALL=100

class Verbose(object):
    """Verbose print in the package
    
    Attributes:
        VERBOSITY: int, default = 0:
            Level of verbosity.
            
            Verbosity levels:
                SIMPLE: Simple messages.
                SYSTEM: System operations.
                VERIFY: Message to verify operations
                DEEP: Deep debugging messages
                ALL: All debugging messages
                
    Methods:
        print(level,msg):
            Print a message if level<=VERBOSITY.
    
    Example:
    
        Verbose.print(1,"Hello world") #No print
        
        Verbose.print(0,"Hello world") #Print

        Verbose.VERBOSITY=1
        Verbose.print(1,"Hello world") #Print
        
        Verbose.VERBOSITY=2
        Verbose.print(1,"Hello world") #Print
        
        Verbose.VERBOSITY=2
        Verbose.print(4,"Hello world") #No print
    """
    VERBOSITY=VERB_ALL
    def print(level,*args):
        if level<=Verbose.VERBOSITY:
            print("  "*level+f"VERB{level}::{inspect.stack()[1][3]}::",*args)
            
#Alias
verbose=Verbose.print

class PrynglesCommon(object):
    """Base class of the package.
    
    All major classes are children of PrynglesCommon class.
    """
    def __init__(self):
        pass
    
    def save_to(self,filename):
        """Save object to a binary file
        
        Parameters:
            filename: string:
                Name of the file where the object will be stored.
        
        Notes:
            Based on https://betterprogramming.pub/load-fast-load-big-with-compressed-pickles-5f311584507e.
        """
        verbose(VERB_SYSTEM,f"Saving object to {filename}")
        pikd = open(filename,"wb")
        dill.dump(self, pikd)
        pikd.close()
            
    def load_from(self,filename):
        """Read object from a binary file.
        
        Parameters:
            filename: string:
                Name of the file where the object is stored.        
        """
        verbose(VERB_SYSTEM,f"Loading object from {filename}")
        pikd = open(filename,"rb")
        data = dill.load(pikd)
        pikd.close()
        verbose(VERB_VERIFY,f"Transferring data to new object")
        self.__dict__=data.__dict__
        return data
    
    def __str__(self):
        """Show content of an object
        
        This method determines the default behavior of the command:
        
            print(object)
        """
        #Remove private attributes
        return str({k:v for k,v in self.__dict__.items() if k[0]!='_'})

from pryngles.version import *

#Constants
from pryngles.consts import *

#Utility modules
from pryngles.misc import *
from pryngles.extensions import *
from pryngles.science import *
from pryngles.plot import *
from pryngles.orbit import *
from pryngles.scatterer import *

#Legacy module
from pryngles.legacy import *

#Core modules
from pryngles.sampler import *
from pryngles.spangler import *
from pryngles.body import *
from pryngles.system import *

#Reset verbosity
Verbose.VERBOSITY=VERB_NONE

#This aliases does not work in modules
print_df=Misc.print_df
sci=Science
