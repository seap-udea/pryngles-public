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

# # Pryngles module: body 

from pryngles import *

# ## External modules

# ## The body class
# 
# The Body class is one of the most important classes in the package. 

Body_doc="""A general body.  This calss is not intended to be used independently, just for inheritance purposes.
    
Initialization attributes:

    kind : string
        One of the kind of bodies defined in the package (see _BODY_KINDS)

    primary: Body
        Object in the center of the orbit of this body.

    orbit: Props
        Object with the orbital properties of the body (eg. orbit.m is the mass)
        see each specific Body definition for attributes.
        orbit must be compatible with rebound.

    physics: Props
        Object with the physical properties of the body (eg. physics.radius)
        see each specific Body definition for attributes.

    optics: Props
        Object with the optical properties of the body (eg. physics.lamb_albedo)
        see each specific Body definition for attributes.

Secondary attributes:

    hash: string
        Hash of the object, ie. a unique string identifying the object 
        (see hash Python function)

    parent: Body
        Body around this body is.  Attribute parent is the same as body.

    childs: list
        List with child bodies (bodies which have this body) as the center.

Public methods:

    update_body(props):
        Update a given property.
"""

class BodyDefaults(object):
    """
    These are the default attributes for any body.
    
    orbit:
    
        These attributes should be compatible with rebound.
    
        m: float [rebound mass units], default = 1
            Mass of the body.  If m = 0 the body does not produce gravitation.
            
    physics:
    
        radius: float [rebound length units], default = 1
            Radius of the body.
            
    optics:
            
        nspangles: int, default = 1000
            Number of spangles on which the object will be discretized.
    """
    orbit=dict(m=1)
    physics=dict(radius=1)
    optics=dict(nspangles=1000)

BODY_KINDS=[]
class Body(PrynglesCommon):
    
    def __init__(self,defaults,kind,primary,orbit,physics,optics):
        
        #Update default properties
        new_orbit=defaults.orbit.copy()
        new_physics=defaults.physics.copy()
        new_optics=defaults.optics.copy()
        
        new_orbit.update(**orbit)
        new_physics.update(**physics)
        new_optics.update(**optics)
        
        self.kind=kind
        self.primary=primary
        if self.primary is not None:
            if not isinstance(self.primary,Body):
                raise AssertionError(f"Primary is not a valid Object")
            else:
                self.parent=primary
                primary._update_childs(self)
        
        self.orbit=Props(**new_orbit)
        self.physics=Props(**new_physics)
        self.optics=Props(**new_optics)
        if 'hash' in self.orbit.__dict__:
            self.hash=self.orbit.hash
        else:
            self.hash=str(hash(self))
    
    def update_body(self,**props):
        """Update properties of the Body.
        
        Parametes:
            orbit: Props                
            physics: Props
            optics: Props
                Properties to update. The current object orbit is updated with new 
                values provided in this new object
                
        Example:
        
            B.update_body(orbit=dict(m=2))
                This only update the attribute m of orbit.
        """
        for prop in PROP_TYPES:
            if prop in props and type(props[prop]) is dict:
                new_prop=self.__dict__[prop].__dict__
                new_prop.update(**props[prop])
                props[prop]=Props(**new_prop)
        
        self.__dict__.update(props)
        self._update_childs()
        self._update_parent()
    
    def _update_childs(self,child=None):
        if 'childs' not in self.__dict__:
            self.childs=[]
        if child is not None:
            self.childs+=[child]
            
    def _update_parent(self,parent=None):
        if 'parent' not in self.__dict__:
            self.parent=parent
        elif parent is not None:
            self.parent=parent
            parent._update_childs(self)

Body.__doc__=Body_doc

# ## Testing


