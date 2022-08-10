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
import unittest
from pryngles import *
class Test(unittest.TestCase):
    def test_common(self):
        
        import numpy as np
        import pandas as pd
        import os
        
        p=PrynglesCommon()
        p.casa=dict(perro=0,gato=3)
        p.data=pd.DataFrame(np.random.rand(4000,100))
        p.save_to("/tmp/save.pck")
        print("File size:",os.path.getsize("/tmp/save.pck")/1e6," Mb")
        g=PrynglesCommon()
        g.load_from("/tmp/save.pck")
        print(g.casa,np.array(g.data).shape)
        
        p.save_to("/tmp/save.bz2",compressed=True)
        print("File size:",os.path.getsize("/tmp/save.bz2")/1e6," Mb")
        g=PrynglesCommon()
        g.load_from("/tmp/save.bz2",compressed=True)
        print(g.casa,np.array(g.data).shape)
                
    def test_fun(self):
        
        #Get path
        filepath=Misc.get_data("diffuse_reflection_function.data")
        print(filepath)
        
        #print_df dataframe
        import pandas as pd
        import numpy as np
        df=pd.DataFrame(np.zeros((5,3)),columns=["a","b","c"])
        Misc.print_df(df)
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
