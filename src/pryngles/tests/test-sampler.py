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
    def test_circle(self):

        #Generate circle
        S=Sampler(1000,seed=10)
        S.gen_circle()
        S.plot()
        S.plot(c='b',spangled=dict(color='r'))
        S.ax.set_title(f"N = {S.N}, dmed = {S.dmed:.4f}, deff = {S.deff:.4f}",fontsize=10)
        S.fig.tight_layout()
        
    def test_ring(self):

        #Generate rings
        S=Sampler(500,seed=10)
        S.gen_ring([[0.0,0.3],[0.5,0.1]])
        print(S.aes)
            
        #Test area
        print(S.A)
        print(np.pi*(1-0.3**2)-np.pi*(0.6**2-0.5**2))
    
        S.plot()
        S.plot(spangled=dict(color='r'))
        S.ax.set_title(f"N = {S.N}, dmed = {S.dmed:.4f}, deff = {S.deff:.4f}",fontsize=10)
        S.fig.tight_layout()

    def test_sphere(self):
        
        #Generate sphere
        S=Sampler(100,seed=10)
        S.gen_sphere()
        S.plot()
        S.plot(spangled=dict(color='r'))
        S.ax.set_title(f"N = {S.N}, dmed = {S.dmed:.4f}, deff = {S.deff:.4f}",fontsize=10)
        S.fig.tight_layout()
        
    def test_purge(self):
        
        #Generate sphere
        S=Sampler(1000,seed=10)
        S.gen_sphere()
        S.purge_sample()
        S.plot()
        S.ax.set_title(f"N = {S.N}, dmed = {S.dmed:.4f}, deff = {S.deff:.4f}",fontsize=10)
        S.fig.tight_layout()
        
    def test_pre(self):
        sp=Sampler(preset="sphere",N=2750)
        print(sp.Npreset,sp.N)
        sp.plot(spangled=dict(color='b',alpha=0.1))
    

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
