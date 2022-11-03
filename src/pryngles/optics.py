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
#!/usr/bin/env python
# coding: utf-8

# # Pryngles module: optics
# 
# Template of a module

# ## External modules

#@external
from pryngles import *
#@end:external

Scatterer_doc="""This is the basic class of a scatterer
"""

# ## Class: Scatterer

#General modules
from pryngles.plot import *
from pryngles.orbit import *

#@class
class Scatterer(PrynglesCommon):
    
    def __init__(self,system):
        self.sys = system
        self.test = True
        
    def compute_angles(self):
        """
        """
        for bhash,body in self.sys.bodies.items():
            if body.kind == "Planet":
                center = body.center_ecl
                        
        azim,incli = Science.spherical(self.sys.n_obs)[1:]
        Rx = self.rotation_matrix_x(np.pi/2-incli)
        Rz = self.rotation_matrix_z(np.pi/2-azim)
        
        luz_equ = spy.unorm(center-self.sys.center_root)[0]
        luz_obs = np.matmul(Rx, np.matmul(Rz, luz_equ))
        
        self.phase_angle = np.dot(luz_obs,np.array([0,0,1]))
        
        for bhash,body in self.sys.bodies.items():
            
            if body.kind == "Star":
                verbose(VERB_SIMPLE,f"Body is a star... skipping")
                continue
                
            elif body.kind == "Planet":
                self.etaps = body.sg.data.cos_luz
                self.zetaps = body.sg.data.cos_obs
                t1 = self.phase_angle - self.zetaps*self.etaps
                t2 = np.sin(np.arccos(self.etaps))*np.sin(np.arccos(self.zetaps))
                t3 = t1/t2
                t3[t3 > 1] = 1.0
                t3[t3 < -1] = -1.0
                self.phidiffps = np.pi - np.arccos(t3)
                self.phidiffps[abs(t2) <= 1e-9] = 0.0 
                self.phidiffps[body.sg.data.y_obs < 0] *= -1
                
            elif body.kind == "Ring":
                self.etars = body.sg.data.cos_luz
                self.zetars = body.sg.data.cos_obs
                t1 = self.phase_angle - self.zetars*self.etars
                t2 = np.sin(np.arccos(self.etars))*np.sin(np.arccos(self.zetars))
                t3 = t1/t2
                t3[t3 > 1] = 1.0
                t3[t3 < -1] = -1.0
                self.phidiffrs = np.pi - np.arccos(t3)
                self.phidiffrs[abs(t2) <= 1e-9] = 0.0 
                
            else: 
                continue

    
    def rotation_matrix_x(self,angle):
        Rm = np.array([[1,0,0],[0,np.cos(angle), np.sin(angle)],[0,-np.sin(angle),np.cos(angle)]])
        return Rm
    
    def rotation_matrix_z(self,angle):
        Rm = np.array([[np.cos(angle), -np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
        return Rm
        
#@end:class        

Scatterer.__doc__=Scatterer_doc


# ### The end

