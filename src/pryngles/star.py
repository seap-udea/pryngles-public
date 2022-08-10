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

# # Pryngles module: stars

from pryngles import *

# ## External modules

sci=Science
print_df=Misc.print_df

# ## Star Class

class StarDefaults(object):
    """
    These are the default attributes for bodies of the kind 'Star'.
    
    DEVELOPER:
        You may add new attributes as the model gets more complex.
        Please document properly each attribute.

    """
    orbit=dict(
        m=1,
    )
    physics=dict(
        radius=1,
        prot=1,
        i=0,#Inclination of the rotational axis
        roll=0,
        alpha=0,#Zero meridian
        t0=0,
    )
    optics=dict(
        nspangles=1000,
        limb_coeffs=[],
    )

BODY_KINDS+=["Star"]
class Star(Body):
    """A star.

    Initialization attributes:
        
        primary: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.
            If None the object is the center of the orbit specification for other objects.
        
        orbit:

            These attributes should be compatible with rebound.

            m: float [um], default = 1: 
                Mass of the star. It should be different than zero.

        physics:

            radius: float [ul], default = 1:
                Radius of the star.

            prot: float [ut], default = 1:
                Period of rotation of the star.
                
            i: float [rad], default = 0:
                Inclination of the ring with respect to the ecliptic plane.

            roll: float [rad], default = 0:
                Roll angle.  This is the angle with respect to ecliptic x-axis in which 
                the normal to the ring plane is rotated.
                
            alpha_equ: float [rad], default = 0:
                Longitude of the zero meridian of the object.
                
            t0: float [ut], default = 0:
                Initial time for zero meridian.

        optics:

            limb_coeffs: list [adimensional], default = []
                List of limb darkening fit coefficients.  See Science.calc_limbdarkening.
                
                Models in: https://pages.jh.edu/~dsing3/David_Sing/Limb_Darkening.html
                Coefficients available at: https://pages.jh.edu/~dsing3/LDfiles/LDCs.CoRot.Table1.txt

            nspangles: int, default = 1000:

                Number of spangles on which the star will be discretized.

    Derived attributes:
    
        physics:
            
            wrot: float [rad/ut]:
                Rotational angular velocity.
                
            n_equ: array(3):
                Rotational axis vector.
    
    Methods:
    
        update_body(**pars):

            This method compute some derived attributes like.

    Notes:

        See Body class documentation.
    
    """
    def __init__(self,
                 primary=None,
                 orbit=StarDefaults.orbit,
                 physics=StarDefaults.physics,
                 optics=StarDefaults.optics
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,StarDefaults,"Star",primary,orbit,physics,optics)

        #Check primary
        if self.primary is not None:
            if self.primary.kind=="Planet":
                raise ValueError(f"Planet cannot be the primary of a Star")
                
        #Update properties
        self.update_body(**self.__dict__)
        
    def update_body(self,**pars):
        Body.update_body(self,**pars)
        
        #Update physics
        
        #Rotational angular velocity
        self.physics.wrot=2*np.pi/self.physics.prot
        
        #Rotation axis:
        self.physics.n_equ,one=spy.unorm(
            sci.cartesian([1,self.physics.roll,90*Consts.deg-self.physics.i])
        )


def spangle_body(self,seed=0):
    """
    Spangle the surface of the star
    """
    
    #Create spangler
    self.sp=Spangler(
        nspangles=self.optics.nspangles,
        body_hash=self.hash,
        spangle_type=STAR_SPANGLE,
        n_equ=self.physics.n_equ,
        alpha_equ=self.physics.alpha,
        w_equ=self.physics.wrot,
        t0_equ=self.physics.t0,
    )
    
    #Populate spangler
    self.sp.populate_spangler(
        scale=self.physics.radius,
        seed=seed,
        geometry="sphere",        
    )

Star.spangle_body=spangle_body


