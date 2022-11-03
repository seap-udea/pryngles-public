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
	def test_interface(self):
	    
	    Verbose.VERBOSITY=VERB_NONE
	    
	    class MySurface(Scatterer):
	        def __init__(self,phase_law,**params):
	            self.A=params["A"]
	
	        def get_albedo(self,eta,zeta,delta,lamb,**params):
	            albedo=self.A*eta
	            return albedo
	        
	    S=MySurface(None,A=0.5)
	    print(S.get_albedo(0.5,0,0,0))
	    
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_lambsurface(self):
	    
	    global LA
	    
	    Verbose.VERBOSITY=VERB_NONE
	    
	    LA=LambertianGraySurface(AL=0.5)
	    etas=np.linspace(0,1,1000)
	    fig,axs=plt.subplots(1,1)
	
	    ax=axs
	    ax.plot(etas,LA.get_albedo(etas,0,0,0))
	    ax.set_xlabel(r"$\zeta = \cos Z$")
	    ax.set_ylabel(r"$\alpha$")
	    ax.set_title(rf"Planetary Lambertian Albedo, $A_L=${LA.AL}");
	
	    fig.tight_layout()
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_lambatmos(self):
	    
	    global LA
	    
	    Verbose.VERBOSITY=VERB_NONE
	    
	    LA=LambertianGrayAtmosphere(AS=0.5)
	    etas=np.linspace(0,1,1000)
	    fig,axs=plt.subplots(1,1)
	
	    ax=axs
	    ax.plot(etas,LA.get_albedo(etas,0,0,0))
	    ax.set_xlabel(r"$\zeta = \cos Z$")
	    ax.set_ylabel(r"$\alpha$")
	    ax.set_title(rf"Atmospheric Lambertian Albedo, $A_S=${LA.AS}");
	
	    fig.tight_layout()
	    Verbose.VERBOSITY=VERB_NONE
	
if __name__=="__main__":
    unittest.main(argv=['first-arg-is-ignored'],exit=False)
    