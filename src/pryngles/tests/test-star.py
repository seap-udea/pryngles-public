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
    def test_star(self):
        S=Star()
        print(S.physics)
        print(S.hash)
        
        #Check derived properties
        self.assertEqual(np.isclose([S.physics.wrot],
                                    [2*np.pi/StarDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        
        S.update_body(physics=dict(m=2))
        print(S.physics)
        
        #Check exception: primary could not be different from None or Body
        self.assertRaises(AssertionError,lambda:Star(primary="Nada"))
        
    def test_sp(self):
        S=Star()
        S.spangle_body()
        print_df(S.sp.data.tail())
        S.sp.plot3d()
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
