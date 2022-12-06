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
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# External required packages
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

from pryngles import *

import spiceypy as spy
import numpy as np
from copy import deepcopy
from anytree import NodeMixin,RenderTree


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Body
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Body(Orbody):
    """A general body.  This class is not intended to be used independently, just for inheritance purposes.
        
    Initialization attributes:
    
        kind : string:
            One of the kind of bodies defined in the package (see _BODY_KINDS)
            Defined objects are: "Star", "Planet", "Ring".
    
        defaults : OrderedDict:
            Dictionary with the properties of the object.
    
        parent: Class Body:
            Object in the center of the orbit of this body.
    
        **properties: dicitionary:
            Specification of the body properties.  All objects of the class Body has the following
            properties by default:
            
            name: string, default = None:
                Name of the object, ie. a unique string identifying the object.  It can be provided
                by the user or automatically set by the initializer using a unique hash 
                (see hash Python function).
    
            orbital properties: 
                Object with the orbital properties of the body (eg. orbit.m is the mass)
                see each specific Body definition for attributes.
                orbit must be compatible with rebound.
    
                    m: float [rebound mass units], default = 1:
                        Mass of the body.  If m = 0 the body does not produce gravitation.
    
            physical properties:
    
                Object with the physical properties of the body (eg. physics.radius)
                see each specific Body definition for attributes.
    
                    radius: float [rebound length units], default = 1:
                        Radius of the body.
    
                    prot: float [ut], default = 1:
                        Period of rotation of the star.
    
                    i: float [rad], default = 0:
                        Inclination of the body equator with respect to the ecliptic plane.
    
                    roll: float [rad], default = 0:
                        Roll angle.  This is the angle with respect to ecliptic x-axis in which 
                        the normal to the object equatorial plane is rotated.
    
                    alpha_equ: float [rad], default = 0:
                        Longitude of the zero meridian of the object.
    
                    q0: float [ut], default = 0:
                        Initial longitude for zero meridian.
    
            optical properties:
    
                Object with the optical properties of the body (eg. physics.lamb_albedo)
                see each specific Body definition for attributes.
    
                    nspangles: int, default = 1000:
                        Number of spangles on which the object will be discretized.
                        
                    spangle_type: int, default = SOLID_SPANGLE:
                        Type of spangles of the body.
                        
                    preset: boolean, default = True:
                        If True spangle object from a preset.
    
    Derived attributes:
    
            wrot: float [rad/ut]:
                Rotational angular velocity.
    
            n_equ: array(3):
                Rotational axis vector in the ecliptic system.
        
    Secondary attributes:
    
        childs: list
            List with child bodies (bodies which have this body) as the center.
    
    Public methods:
    
        update_body(**props):
            Update a given set of properties.
            
    Examples:
    
        Create a body with None parent and name = 'B':
        
            B=Body("Body",BODY_DEFAULTS,None,name='B',m=2,c=2)
            
        Create a body having parent the Body "B" defined before:
             
            C=Body("Body",BODY_DEFAULTS,B,name="C")
    
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    def __init__(self,kind,defaults,parent,**props):

        #Kind, parent and child attributes
        self.kind=kind
        self.__defaults=defaults
        
        #Prepare key attributes
        self.sg=None

        #Name of the object
        if 'name' in props:
            name=self.name=str(props["name"])
        elif 'name_by_kind' in props:
            name=self.name=self.kind
        else:
            name=self.name=str(hash(self))

        #Legacy
        if 'primary' in props:
            parent=props["primary"]
        if 'optics' in props:
            props.update(props["optics"])
        if 'orbit' in props:
            props.update(props["orbit"])
        if 'physics' in props:
            props.update(props["physics"])

        #Update childs and parent
        if parent is not None:
            if not isinstance(parent,Body):
                raise AssertionError(f"Parent is not a valid Object: {type(parent)}, {isinstance(parent,Body)}")
            else:
                self.parent=parent
                parent._update_childs(self)

        #Update parent and childs        
        self._update_parent(parent)
        self._update_childs()

        #Update default properties
        self.__dict__.update(defaults)
        #Set name
        self.name=name
        #Update body
        self.update_body(**props)
    
    def update_body(self,**props):
        """Update properties of the Body.
        
        Parametes:
            **props: dictionary:
                Properties to update. The current object is updated with new 
                values provided in this new object
                
        Example:
            B.update_body(m=2)
                This only update the attribute m of orbit.
        """
        for prop in props:
            if prop in self.__defaults or prop in REBOUND_ORBITAL_PROPERTIES:
                self.__dict__[prop]=props[prop]
            else:
                raise ValueError(f"Property {prop} not identified in object {self.kind}")
                
        self.elements={k:v for k,v in self.__dict__.items() if k in REBOUND_ORBITAL_PROPERTIES}
        
        verbose(VERB_VERIFY,"Updating Body")
        self._update_properties()
    
    def _update_childs(self,child=None):
        if 'childs' not in self.__dict__:
            self.childs=dict()
        if child is not None:
            verbose(VERB_VERIFY,f"Add child {child.name} to body {self.kind} ({self.name})")
            self.childs[child.name]=child
            
    def _update_parent(self,parent=None):
        if 'parent' not in self.__dict__:
            if parent:
                verbose(VERB_VERIFY,f"Add parent {parent.name} to body {self.kind} ({self.name})")
            self.parent=parent
        elif parent is not None:
            verbose(VERB_VERIFY,f"Add parent {parent.name} to body {self.kind} ({self.name})")
            self.parent=parent
            parent._update_childs(self)
    
    def _update_properties(self):
        verbose(VERB_VERIFY,"Updating properties of Body")
        #Rotational angular velocity
        self.wrot=2*np.pi/self.prot
        #Rotation axis
        self.n_equ=sci.cartesian([1,self.roll,90*Consts.deg-self.i])
    
    def show_tree(self):
        print(RenderTree(self))
        

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Tested methods from module file body
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def spangle_body(self):
        """
        Spangle the surface of the body
        """
        
        #Create spangler
        self.sg=Spangler(
            nspangles=self.nspangles,
            name=self.name,
            n_equ=self.n_equ,
            alpha_equ=self.alpha,
            w=self.wrot,
            q0=self.q0,
        )
        
        #Populate spangler
        self.sg.populate_spangler(
            shape=self.shape,
            spangle_type=self.spangle_type,
            scale=self.radius,
            seed=self.seed,
            preset=self.preset,
            **self.geometry_args,
        )
        
        #Additional properties in the Spangler DataFrame
        if self.kind=="Star":
            self.sg.data.source=True
        
        self.sg.set_observer()
        self.sg.set_luz()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Star
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Star(Body):
    """A star.

    Initialization attributes:
        
        parent: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.

            If None the object is the center of the orbit specification for other objects.
            
            Object parent for a star should be another star.
        
        **props: dictionary:
            List of properties for star.  For the complete set of default values of the properties
            see STAR_DEFAULTS.  Description of properties are available in the Body class documentation.
            
            Additional properties:
            
                limb_coeffs: list [adimensional], default = []:
                    List of limb darkening fit coefficients.  See Science.calc_limbdarkening.

                    Models in: https://pages.jh.edu/~dsing3/David_Sing/Limb_Darkening.html
                    Coefficients available at: https://pages.jh.edu/~dsing3/LDfiles/LDCs.CoRot.Table1.txt
                    
                spangle_type: int, default = STAR_SPANGLE:
                    Type of spangles

    Derived attributes:
    
    Methods:
    
        update_body(**pars):

            This method compute some derived attributes like.

    Notes:

        See Body class documentation.
    
    """
    def __init__(self,
                 parent=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Star",STAR_DEFAULTS,parent,**props)

        #Check parent
        if self.parent is not None:
            if self.parent.kind!="Star":
                raise ValueError(f"Only another Star can be the parent of a Star (you provided {self.parent.kind})")

        self._update_star_properties()
        
    def _update_star_properties(self):
        """Update specific properties of the star
        
        Properties to update:
        
            norm_limb_darkening: float:
                Limb darkening function normalization.
                Requires: limb_coefs.

        """
        verbose(VERB_VERIFY,"Updating properties of Star")
        
        #Compute limbdarkening at r = 0 to initialize normalization constant
        sci.limb_darkening(0,self.limb_coeffs)
        
        #Store limb darkening normalization
        self.norm_limb_darkening=SCIENCE_LIMB_NORMALIZATIONS[hash(tuple(self.limb_coeffs))]
        
    def update_star(self,**props):
        """General update propeties of the Star
        """
        verbose(VERB_VERIFY,"Updating star")
        
        Body.update_body(self,**props)
        self._update_star_properties()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Planet
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Planet(Body):
    """A planet.

    Initialization attributes:
        
        parent: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.
            If None the object is the center of the orbit specification for other objects.

        **props: dictionary:
            List of properties for star.  For the complete set of default values of the properties
            see STAR_DEFAULTS.  Description of properties are available in the Body class documentation.
            
            Additional properties:
            
                x,y,z: float [ul], default = 1.0, 0.0, 0.0:
                    Initial position of the body.

                vy: float, default = 0.0, 1.0, 0.0:
                    Intitial velocity of the body
        
    Derived attributes:
        None.
    
    Notes:

        See Body class documentation.
    
    """
    
    def __init__(self,
                 parent=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Planet",PLANET_DEFAULTS,parent,**props)
        
        #Check parent
        if self.parent is None:
            raise ValueError(f"Parent not provided and it is mandatory for {self.kind}.")
        
        #Update properties
        self.update_planet(**props)

    def _update_planet_properties(self):
        """Update specific properties of the star
        
        Properties to update:
        
            norm_limb_darkening: float:
                Limb darkening function normalization.
                Requires: limb_coefs.

        """
        verbose(VERB_VERIFY,"Updating Planet properties")
        
    def update_planet(self,**pars):
        verbose(VERB_VERIFY,"Updating Planet")
        Body.update_body(self,**pars)
        self._update_planet_properties()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Ring
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Ring(Body):
    """Class Ring.
    
Initialization attributes:
        
        parent: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.
            If None the object is the center of the orbit specification for other objects.

        **props: dictionary:
            List of properties for star.  For the complete set of default values of the properties
            see STAR_DEFAULTS.  Description of properties are available in the Body class documentation.
            
            Additional properties:

            fi: float [adimensional], default = 1:
                Fraction of the radius of the parent object where ring stars.

            fe: float [adimensional], default = 1:
                Fraction of the radius of the parent object where ring ends.

            albedo_gray_normal: float. default = 1: 
                Lambertian (normal) gray (wavelength indpendent) albedo of the spangle.

            tau_gray_optical: float. default = 0:
                Gray (wavelength indpendent) Optical depth of the spangle.  
                If 0 the spangle is entirely opaque to all wavelength, despite its type.            

    Derived attributes:
    
        ri: float:
            Radius of the inner border of the ring in units of the parent radius.

        re: float:
            Radius of the outer border of the ring in units of the parent radius.
            
    Notes:

        See Body class documentation.
    """
    def __init__(self,
                 parent=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Ring",RING_DEFAULTS,parent,**props)
        
        #Check parent
        if self.parent is None:
            raise ValueError(f"Parent not provided and mandatory for {self.kind}.")
        
        #Update properties
        self.update_ring(**props)

    def _update_ring_properties(self):
        """Update specific properties of the star
        
        Properties to update:
        
            ri, re: float:
                Radius of the inner (outer) border of the ring in units of the parent radius.
                Requires: limb_coefs.
                
            radius: float:
                Object radius.
                
            geometry_args: dictionary:
                
        """
        verbose(VERB_VERIFY,"Updating Ring properties")
    
        #Update radius
        self.ri=self.fi*self.parent.radius
        self.re=self.fe*self.parent.radius
        self.radius=self.re
        
        #Update geometry args for spangling purposes
        self.geometry_args=dict(ri=self.ri/self.re)
        
    def update_ring(self,**pars):
        verbose(VERB_VERIFY,"Updating Ring")
        Body.update_body(self,**pars)
        self._update_ring_properties()   



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Observer
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Observer(Body):
    """This class is intended only for legacy purposes.
    """
    def __init__(self,
                 parent=None,
                 **props
                ):
        Body.__init__(self,"Observer",OBSERVER_DEFAULTS,parent,**props)
