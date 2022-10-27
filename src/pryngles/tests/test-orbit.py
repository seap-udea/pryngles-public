##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#

# PlanetaRY spanGLES                                             #
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################
# Main contributors:                                             #
#   Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado         #
##################################################################
import unittest
from pryngles import *
class Test(unittest.TestCase):
	def test_orbody(self):
	
	    S=Orbody(m=1)
	    P=Orbody(parent=S,m=0.1,a=1,e=0.5)
	
	    print(S.m,S.elements)
	    print(P.m,P.elements)
	
	    S.show_tree()
	
	def test_Orbit(self):
	
	    global S4
	
	    #Quintuple system
	    S1=Orbit(m1=1,m2=1,a=1,e=0.7,M=0)
	    S2=Orbit(m1=1,m2=1,a=1,e=0,M=0)
	    S3=Orbit(m1=S1,m2=S2,a=5,e=0)
	    S4=Orbit(m1=S3,m2=1,a=20,e=0,E=45*Consts.deg)
	    S4.calculate_orbit()
	    print(S4.get_states())
	    print(S4.Ps)
	
	    #Using custom units
	    #Quintuple system
	
	    #Initialize positions
	    units=['au','msun','yr']
	    hn=Orbit(
	        m1=1,
	        m2=Orbit(m1=1e-3,m2=1e-7,a=0.5,e=0.0,units=units),
	        units=units,
	        a=20,e=0.0)
	    hn.calculate_orbit()
	    sim,states=hn.get_states()
	    print(states)
	
	    #SImple 
	
	    #Use this code to animate:
	    #Plot.animate_rebound(S4.sim)
	
if __name__=="__main__":
    unittest.main(argv=['first-arg-is-ignored'],exit=False)
    