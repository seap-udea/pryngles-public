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

from pryngles import *
from collections import OrderedDict as odict
from collections.abc import Iterable
import inspect
import os
from sys import maxsize as HASH_MAXSIZE


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Misc
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Misc(object):
    """
    Miscelaneous routines.
    
    This is a set of util routines intended for a diversity of purposes.
    
    Routines included:
    
        get_data(file)
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        """Print DataFrame.
        
        Parameters:
            df: Pandas DataFrame:
                DataFrame to print.
        """
        display(HTML(df.to_html()))
        
    def flatten(collection):
        """Flatten a list of objects

        Examples:
            list(Misc.flatten(["cosa"]))
            list(Misc.flatten([["cosa"]]))
            list(Misc.flatten([["cosa","perro"]]))
            list(Misc.flatten([[1,"perro"],object,float]))
        """
        for i in collection:
            if isinstance(i, Iterable) and not isinstance(i, basestring):
                for subc in Misc.flatten(i):
                    yield subc
            else:
                yield i
                
    def get_methods(my_class):
        """Get a list of the methods for class my_class
        """
        return sorted([member[0] for member in inspect.getmembers(my_class) if '__' not in member[0]])
    
    def calc_hash(obj):
        if type(obj) is dict:
            hash_obj=frozenset(obj.items())
        else:
            hash_obj=obj
        hash_val=str(hash(hash_obj)%((HASH_MAXSIZE+1)*2))
        return hash_val

