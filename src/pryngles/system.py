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

# # Pryngles module: System

from pryngles import *

# ## External modules

import rebound as rb

# ## System Class
# 
# This is the most important class in the whole package.  This class allows to create the planetary system and manipulate it.

System_doc=f"""
Creates a planetary system.

Examples:

    sys=System()
    sys=System(units=['km','msun','s'])
    S=Star()
    P=Planet(primary=S)
    sys=System(stars=S,planets=P)

Initialization attributes:

    units: list of strings, default = ['au','msun','yr']:
        Units used in calculations following the conventions and signs of rebound.
        The order SHOULD always be: length, mass, time.
        
    stars: list or single Body.  default = None: 
        Star(s) in the system.
        
    planets: list or single Body.  default = None: 
        Planet(s) in the system.
        
    rings: list or single Body.  default = None: 
        Ring(s) in the system.

    observers: list or single Body.  default = None: 
        Observer(s) in the system.
    
    rebound: bool, default = True:
        Set True if you want to simulte the orbit of objects (for instance if you want to calculate
        TTVs or TDVs) or False if you want to calculate orbits using Keplerian dynamics (no TTVs or TDVs).
        
        NOTE: Current version ({version}) does not implement yet rebound simulations

Secondary attributes:

    hashes: list
        List of hashes of bodies in the system.
        
    stars, planets, rings, observers: list
        List of the corresponding kind of object in the system.
        
    nstars, nplanets, nrings, nobservers: int
        Number of each kind of body in the system.
    
    nbodies: int
        Number of bodies in the system.
        
    ul, um, ut: float [SI units]
        Value of the conversion factors for each unit.
        
    G: float [ul^3/ut^2/um]
        Value of the gravitational constant.

Usefule private attributes:

    _sim: rebound Simulation.
        Rebound simulation object.
""";

class System(PrynglesCommon):
    
    def __init__(self,
                 units=['au','msun','yr'],
                 stars=None,planets=None,
                 rings=None,observers=None,
                 rebound=False
                ):
        
        #Behavior
        self.rebound=rebound
        self.rebound=False # Rebound not implemented yet
        self._update_rebound(units)

        #Initialize rebound
        self.update_units(units)
        
        #Initialize list of components
        self.hashes=dict()
        for kind in BODY_KINDS:
            lkind=kind.lower()
            exec(f"self.{lkind}s=[]")
            exec(f"self._update_objects('{lkind}s','{kind}',{lkind}s)")

        #Update system
        self._update_system()
    
    def _update_objects(self,attr,kind,comps):
        """Update the list of objects
        """

        if comps is not None:
            #Check if comps is a a list
            try:
                comps[0]
            except:
                comps=[comps]
            
            for comp in comps:
                #Check if primary bodies are already in lists
                if (comp.primary is not None) and (comp.primary.hash not in self.hashes):
                    raise AssertionError(f"Primary of {kind} body is not yet in the system.")
                else:
                    exec(f"self.{attr}+=[comp]")
                    if comp.kind!=kind:
                        raise AssertionError(f"You are attempting to add {kind} with a {comp.kind}")
                    self.hashes[comp.hash]=comp

    def _update_system(self):
        """Update the global properties of the system.
        """
        for kind in BODY_KINDS:
            exec(f"self.{kind}s=0")
        
        pdict=dict()
        for obj in BODY_KINDS:
            lobj=obj.lower()
            exec(f"{obj}=len(self.{lobj}s)",locals(),pdict)
            self.__dict__[f"n{lobj}s"]=pdict[obj]
        self.nbodies=self.nstars+self.nplanets+self.nrings
        
    def update_units(self,units):
        """Update units of the system
        """
        self.units=units
        self._ul,self._um,self._ut=self.units
        self._sim.units=self.units
        self.ul=rb.units.convert_length(1,self._ul,"m")
        self.um=rb.units.convert_mass(1,self._um,"kg")
        self.ut=np.sqrt(self._sim.G*self.ul**3/(self.um*GSI))

    def _update_rebound(self,units):
        """Crteate and update rebound simulation object
        """
        #Units initialization
        self._sim=rb.Simulation()
        self._sim.units=units

System.__doc__=System_doc


def add(self,kind=None,primary=None,orbit=None,physics=None,optics=None):
    """Add an object to the system
    
    Examples:
    
        sys=System()
        S=sys.add(kind="Star",orbit=dict(m=2))
    
    Parameters:
    
        kind: string
            Kind of object: Star, Planet, Ring, Observer.
    
        primary: Body
            Primary object.
    
        orbit: dictionary
            Set of orbital properties (see corresponding body documentation).
            
        physics: dictionary
            Set of physical properties (see corresponding body documentation).
        
        optics: dictionary
            Set of optical properties (see corresponding body documentation).
            
    Returns:
        
        Body
            Body added to the system.
    """
    if kind is None:
        raise AssertionError("You must provide a valid object kind (Star, Planet, Ring, Observer).")

    if kind not in BODY_KINDS:
        raise ValueError(f"Object kind '{kind}' is not recognized.")

    #p.e. kind = 'Star', lkind = 'star'
    kind=kind.capitalize()
    lkind=kind.lower()

    """
    #This code generalize the procedure:
    if orbit is None:
        orbit = StarDefaults.orbit.copy()
    if physics is None:
        physics = StarDefaults.physics.copy()
    if optics is None:
        optics = StarDefaults.optics.copy()
    """
    pdict=dict()
    plist=[]
    for prop in PROP_TYPES:
        exec(f"{prop}={prop}",locals(),pdict)
        if pdict[prop] is None:
            exec(f"{prop}={kind}Defaults.{prop}.copy()",globals(),pdict)
        plist+=[pdict[prop]]

    """
    #This code generalize the procedure
    S=Star(primary=primary,orbit=orbit,physics=physics,optics=optics)
    self._list2Objects("stars","Star",S)            
    """
    obj=eval(f"{kind}(primary=primary,orbit=plist[0],physics=plist[1],optics=plist[2])")
    self._update_objects(lkind+"s",kind,obj)
    self._update_system()
    return obj
    
System.add=add


def remove(self,body_hash):
    """Remove a body from a system.

    Example:
        sys=System()
        S=sys.add(kind="Star",orbit=dict(m=2))
        sys.remove(body_hash=S.hash)
        
    Parameters:
        body_hash: string
            Hash of the body to remove
        
    Remove eliminate body and all the childs and the childs of the childs.
    """
    if body_hash in self.hashes:
        obj=self.hashes[body_hash]
        lkind=obj.kind.lower()

        #Remove child objects
        for child in obj.childs:
            if child.hash in self.hashes:
                self.remove(child.hash)

        #Remove object
        exec(f"self.{lkind}s.remove(obj)")

        #Remove hash from list
        self.hashes.pop(body_hash)

        #Update system
        self._update_system()
    else:
        raise ValueError("No object with hash 'body_hash' in the system")
System.remove=remove


