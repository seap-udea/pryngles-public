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
    def test_temp(self):
        
        #Test
        Science.template()
        """
        self.assertEqual(self.P.Nr,8,True)
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        self.assertRaises(AssertionError,lambda:Observer(primary="Nada"))
        """
        
    def test_coords(self):
        
        #Test spherical
        rqf=Science.spherical([1,1,0])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([+1,+1,1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([-1,1,1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([-1,-1,1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([+1,-1,1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([+1,+1,-1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([-1,1,-1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([-1,-1,-1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([+1,-1,-1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)
        
        #Test cartesian
        xyz=Science.cartesian([1,0,0])
        print(xyz) 
        xyz=Science.cartesian([1,45*Consts.deg,45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,135*Consts.deg,45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,225*Consts.deg,45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,315*Consts.deg,45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,45*Consts.deg,-45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,135*Consts.deg,-45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,225*Consts.deg,-45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,315*Consts.deg,-45*Consts.deg])
        print(xyz) 

    def test_rot(self):
        Verbose.VERBOSITY=1
        
        #Test rotation
        Msys2uni,Muni2sys=Science.rotation_matrix([0,0,1],0)
        print(Msys2uni)

        Msys2uni,Muni2sys=Science.rotation_matrix([0,0,-1],0)
        print(Msys2uni)

        Verbose.VERBOSITY=0


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
