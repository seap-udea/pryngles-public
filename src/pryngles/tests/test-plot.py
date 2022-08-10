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
    def test_fun(self):
        fig=plt.figure()
        ax=fig.add_subplot(111,projection='3d')
        ax.plot([0],[0],[0],'ko',markersize=5)
        Plot.circle3d(ax,(0,0,0),0.5,zDir=[1,1,0],fill='None',alpha=0.3)
        #Decoration
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_zlim(-1,1)
        ax.set_box_aspect([1,1,1])
        ax.set_xlabel("$x_{ecl}$")
        ax.set_ylabel("$y_{ecl}$")
        ax.set_zlabel("$z_{ecl}$")
        ax.view_init(elev=15,azim=1)
        Plot.pryngles_mark(ax)
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
