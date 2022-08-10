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

# # Pryngles module

from pryngles import *

# ## External modules

# ## Observer Class

class ObserverDefaults(object):
    """
    These are the default attributes for bodies of the kind 'Observer'.
    
    DEVELOPER:
        You may add new attributes as the model gets more complex.
        Please document properly each attribute.
        
    orbit:
    
        (In current version, Observer body does not have orbit attributes)
        
    physics:
    
        (In current version, Observer body does not have physics attributes)
            
    optics:
    
        distance: float [rebound length units], default = 1
            Distance at which the observer is located
            
        beta: float [radians], default = 0
            Ecliptic latitude.  Elevation with respect to the plane of the planetary 
            system from which the observer is looking at the system.  beta = 0 is for 
            edge-on observations, beta = 90 degrees is for face-on observations.
            
        lamb: float [radians], default = 0
            Ecliptic longitude. Angle with respect to the zero meridian of the planetary 
            system.
            
            NOTE: zero meridian is defined differently depending on the observed object.
                  in the case of a ring, zero meridian is the direction of vernal equinox.
    
        nspangles: int, default = 1000
            Number of spangles on which the object will be discretized.
    """
    orbit=dict()
    
    physics=dict()
    
    optics=dict(distance=1,beta=0,lamb=0)

BODY_KINDS+=["Observer"]

class Observer(Body):
    """Class Observer.
    
    See Body class documentation.
    
    Additional public attributes:
    
        physics.inclination: float [radians]
            Angle of inclination of the normal to the ecliptic in the system with 
            respect to line-of-sight.  When inclination = 0 system is observed face-on.
            When inclination = 90 degress system is observed edge-on.

    Override methods:
    
        update_body(**pars):
            This method compute additional attributes like (see above).
    """
    def __init__(self,
                 primary=None,
                 orbit=ObserverDefaults.orbit,
                 physics=ObserverDefaults.physics,
                 optics=ObserverDefaults.optics
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,ObserverDefaults,"Observer",primary,orbit,physics,optics)
        
        #Check if observer is attached to any object
        """
        self.primary=primary
        if self.primary is not None:
            self.primary._update_childs(self)
            self._update_parent(self.primary)
        """
        
        #Update properties
        self.update_body(**self.__dict__)
        
    def update_body(self,**pars):
        Body.update_body(self,**pars)
        
        #Here place the commands to update this kind of body
        self.optics.inclination = 90*Consts.deg - self.optics.beta


