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
    def test_init(self):
        Verbose.VERBOSITY=1
        print("Basic definition:")
        sg=Spangler(nspangles=3,spangle_type=GRANULAR_SPANGLE)
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Equ->Obs:\n",sg.M_obs2ecl)

        print("\nAnother definition:")
        sg=Spangler(nspangles=3,n_equ=[1,0,0])
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Obs->Ecl:\n",sg.M_obs2ecl)

        print("\nDefinition observer:")
        sg=Spangler(nspangles=3,body_hash="123",n_equ=[0,1,0],n_obs=[1,0,0])
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Obs->Ecl:\n",sg.M_obs2ecl)
        Verbose.VERBOSITY=0

    def test_pop(self):
        Verbose.VERBOSITY=1
        sg=Spangler(nspangles=100)
        sg.populate_spangler(geometry="ring",scale=1,seed=1,gaps=[[0,0.2],[0.5,0.1]],boundary=0)
        print_df(sg.data.head(5))
        sg.set_observer(n_obs=[1,1,1])
        print_df(sg.data.head(5))
        sg=Spangler(nspangles=1000,body_hash="123",n_equ=[1,0,1])
        sg.populate_spangler(geometry="sphere",scale=2,seed=1)
        print(sg.nspangles,sg.sample.N,len(sg.data))
        Verbose.VERBOSITY=0
        
    def test_plot3d(self):
        Verbose.VERBOSITY=0
        plt.close("all")
        sg=Spangler(nspangles=500,body_hash="123",n_equ=[1,1,1])
        sg.populate_spangler(geometry="sphere",scale=2,seed=1)
        sg.set_luz(n_luz=[1,0,0])
        sg.plot3d(spangled=False,factor=1.3,c='y',s=3)
        
        sg=Spangler(nspangles=500,body_hash="123",n_equ=[1,1,1])
        sg.populate_spangler(geometry="ring",scale=2,seed=1,boundary=0)
        sg.set_luz(n_luz=[0,0,-1])
        sg.plot3d(factor=0.4)

        Verbose.VERBOSITY=0

    def test_plotobs(self):
        Verbose.VERBOSITY=0
        sg=Spangler(nspangles=500,body_hash="123",n_equ=[1,1,1],center_ecl=[1,1,1])

        sg.populate_spangler(geometry="sphere",scale=2,seed=1)
        sg.set_observer(n_obs=[1,0,0])
        sg.plot_obs()

        sg.populate_spangler(geometry="circle",scale=2,seed=1,boundary=0)
        sg.plot_obs()

        sg.populate_spangler(geometry="ring",scale=2,seed=1,gaps=[[0,0.2],[0.5,0.1]],boundary=0)
        sg.plot_obs()
        
        Verbose.VERBOSITY=0

        
    def test_join(self):
        Verbose.VERBOSITY=0

        sg1=Spangler(nspangles=1000,body_hash="123",n_equ=[1,0,1])
        sg1.populate_spangler(geometry="ring",scale=2.5,seed=1,gaps=[[0,1.5/2.5]],boundary=0)

        sg2=Spangler(nspangles=1000,body_hash="345",n_equ=[0,0,1])
        sg2.populate_spangler(geometry="sphere",scale=1,seed=1)

        sgj=Spangler(spanglers=[sg1,sg2],n_obs=[1,0,0],n_luz=[-1,-1,-1])

        sgj.plot3d()
        sgj.plot_obs()
        
        Verbose.VERBOSITY=0

    def test_scale(self):
        Verbose.VERBOSITY=0

        sg=Spangler()
        print_df(sg.data)

        sg.set_scale(5)
        print_df(sg.data)

        Verbose.VERBOSITY=0


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
