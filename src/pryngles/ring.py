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
#!/usr/bin/env python
# coding: utf-8

# # Pryngles module: rings

from pryngles import *
sci=Science

# ## External modules

import spiceypy as spy
import math as mh
import numpy as np
import copy

# ## Ring default properties

class RingDefaults(object):
    """
    These are the default attributes for bodies of the kind 'Ring'.
    
    DEVELOPER:
        You may add new attributes as the model gets more complex.
        Please document properly each attribute.
        
    orbit:
    
        (In current version, Ring body does not have orbit attributes)
        
    physics:
    
        fi: float [adimensional], default = 1:
            Fraction of the radius of the primary object where ring stars.
            
        fe: float [adimensional], default = 1:
            Fraction of the radius of the primary object where ring ends.
            
        i: float [radians], default = 0:
            Inclination of the ring with respect to the ecliptic plane.
            
        roll: float [radians], default = 0:
            Roll angle.  This is the angle with respect to ecliptic x-axis in which 
            the normal to the ring plane is rotated.
            
    optics:
    
        nspangles: int, default = 1000: 
            Number of spangles on which the object will be discretized.
            
        albedo_gray_normal: float. default = 1: 
            Lambertian (normal) gray (wavelength indpendent) albedo of the spangle.
            
        tau_gray_optical: float. default = 0:
            Gray (wavelength indpendent) Optical depth of the spangle.  
            If 0 the spangle is entirely opaque to all wavelength, despite its type.            
    """
    orbit=dict()
    
    physics=dict(fi=1.0,fe=2.0,i=0.0,roll=0.0)
    
    optics=dict(nspangles=1000,
                albedo_gray_normal=1,
                tau_gray_optical=0
               )

# ## Ring Class

BODY_KINDS+=["Ring"]

class Ring(Body):
    """Class Planet.
    
    See Body class documentation.
    
    Additional public attributes:
    
        ri: float [rlu]:
            Radius of the inner border of the ring

        re: float [rlu]:
            Radius of the outer border of the ring
            
        spangles: list of Spangle objects:
            List of spangles covering the surface (for spangle options see Spangle documentation)
    
    Override methods:
    
        update_body(**pars):
            This method compute additional attributes like (see above).
    """
    def __init__(self,
                 primary=None,
                 orbit=RingDefaults.orbit,
                 physics=RingDefaults.physics,
                 optics=RingDefaults.optics
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,RingDefaults,"Ring",primary,orbit,physics,optics)
        
        #Check primary
        if self.primary is None:
            raise ValueError(f"Primary not provided and mandatory for {self.kind}.")
        
        #Update properties
        self.update_body(**self.__dict__)


# ### Update body

def update_body(self,**pars):
    Body.update_body(self,**pars)

    #Here place the commands to update this kind of body
    self.ri=self.physics.fi*self.primary.physics.radius
    self.re=self.physics.fe*self.primary.physics.radius

Ring.update_body=update_body


