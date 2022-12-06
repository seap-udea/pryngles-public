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

from abc import ABC, abstractmethod
from scipy.optimize import bisect
from scipy.integrate import quad,dblquad
from scipy.interpolate import interp1d,interp2d


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Scatterer
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Scatterer(PrynglesCommon,ABC):
    """An scatterer surface
    
    Initialization attributes:
    
         params: dictionary:
             Other parameters of the phase law.
    
    Mandatory methods:
    
        __init__(self,phase_law:lambda,**params)->int:
            This method should return a hash of the object.
            
        get_albedo(eta:float,zeta:float,delta:float,lamb:float,**params)->float
            This method must provide the albedo.
        
    Class methods:
    
        register(scatterer,params):
            Register scatterer for future uses.
    
        reset_catalogue():
            Reset the catalogue of scatterers.
            
            Usage: Scatterer.reset_catalogue()
    
    Usage:
        You can create a Scatterer which implements this class:
    
            class MySurface(Scatterer):
                def __init__(self,**params):
                    if self.register(self,params):
                        #Read parameters of the scatterer
                        self.A=params["A"]
                        #Initialize scatterer
                        self._initialize_scatterer()
    
                #Mandatory methods
                def get_albedo(self,eta,zeta,delta,lamb,**params):
                    albedo=self.AA*eta
                    return albedo
    
                # Private methods to prepare scatterer
                def _initialize_scatterer(self):
                    self.AA=self.A**2
    
    
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    @abstractmethod
    def __init__(self,**params)->str:
        pass
    
    @abstractmethod
    def get_albedo(self,eta:float,zeta:float,delta:float,lamb:float,**params)->float:
        pass
    
    @classmethod
    def register(self,scatterer,params):
        """Register scatterer
        """
        scatterer.params=params
        scatterer.params["name"]=scatterer.__class__.__name__
        scatterer.hash=Misc.calc_hash(params)
        if scatterer.hash in SCATTERERS_CATALOGUE:
            verbose(VERB_SIMPLE,f"Scatterer with name {scatterer.params['name']} and hash {scatterer.hash} already exist at {id(SCATTERERS_CATALOGUE)}")
            scatterer.__dict__=deepcopy(SCATTERERS_CATALOGUE[scatterer.hash].__dict__)
            return False
        else:
            verbose(VERB_SIMPLE,f"Creating a new scatterer with name {scatterer.params['name']} and hash {scatterer.hash}")
            scatterer.params["hash"]=scatterer.hash
            SCATTERERS_CATALOGUE[scatterer.hash]=scatterer
            return True
        
    @classmethod
    def reset_catalogue(self):
        """Reset catalogue of scatterers
        """
        SCATTERERS_CATALOGUE=dict()
        


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class NeutralSurface
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class NeutralSurface(Scatterer):
    """Neutral surface.
    """
    def __init__(self,**params):
        if self.register(self,params):
            pass
    
    def get_albedo(self,eta,zeta,delta,lamb,**params):
        return 1
    


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class BlackBodySurface
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class BlackBodySurface(Scatterer):
    """Black body surface
    """
    def __init__(self,**params):
        if self.register(self,params):
            pass
    
    def get_albedo(self,eta,zeta,delta,lamb,**params):
        return 0
    


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class LambertianGraySurface
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class LambertianGraySurface(Scatterer):
    """Lambertian Gray Surface.
    
    This is the scatterer corresponding to a surface having a gray lambertian Albedo.
    
    Parameters:
    
        phase_law: function, default=lambda eta,zeta,delta,lambda:eta :

            Law of reflection (by default is Lambertian, see Russel, 1916)

            The phase_law must obey the following prototype:

                phase_law(eta,zeta,delta,lamb,**params):
                    '''Phase law of the surface

                    Parameters:
                        eta: float:
                            cosine of the incoming angle.

                        zeta: float:
                            cosine of the outgoing angle.

                        delta: float:
                            difference between the incoming and outgoing azimuth.

                        lamb: float:
                            Wavelength.

                        parameters: dictionary: 
                            Other parameters of the phase law.

                    Return:
                        Wavelength dependent albedo.
                    '''
                    ...

                Other law is the Lommel-Seeliger law:

                    phase_law = lambda eta,zeta,delta,params:eta*zeta/(eta+zeta) (see Russel, 1916)

    """
    
    def __init__(self,**params):

        
        if self.register(self,params):
            verbose(VERB_SIMPLE,f"Initializing {self.params['name']} with hash {self.hash}")
            
            #Phase law
            if "phase_law" in params:
                self.phase_law=params["phase_law"]
            else:
                self.phase_law=lambda eta,zeta,delta,lamb,params:eta

            #Gray albedo
            self.AL=params["AL"]

            #Calculate the gammap parameter
            self.gammap0=self._find_gammap()

            #Accelerate the calculation of the albedo
            self._accelerate_lambertian_albedo()

    def get_albedo(self,eta,zeta,delta,lamb,**params):
        return self._get_albedo(eta)
        
    #####################################
    #Complimentary routines
    #####################################
    def _calc_lambertian_albedo(self,eta,gammap0=1):
        if eta==0:return self.AL
        integrand=lambda zeta:self.phase_law(eta,zeta,0,0,0)/eta
        AL=2*np.pi*gammap0*quad(integrand,0,1)[0]
        return AL

    def _find_gammap(self):
        function=lambda gammap0:self._calc_lambertian_albedo(1,gammap0)-self.AL
        gammap0=bisect(function,0.0,1.0,rtol=1e-3)
        return gammap0 if gammap0<=1 else 1
    
    def _accelerate_lambertian_albedo(self):
        etas=np.linspace(0.0,1.0,20)
        ALs=np.array([self._calc_lambertian_albedo(eta,gammap0=self.gammap0) for eta in etas])
        self._get_albedo=interp1d(etas,ALs)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class LambertianGrayAtmosphere
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class LambertianGrayAtmosphere(Scatterer):
    """Lambertian Gray Atmopshere.
    
    This is the scatterer corresponding to plane-parallel gray lambertian atmosphere
    """
    
    def __init__(self,**params):

        if self.register(self,params):        
            verbose(VERB_SIMPLE,f"Initializing {self.params['name']} with hash {self.hash}")
            
            #Gray albedo
            self.AS=params["AS"]

            #Load reflection functions
            self._load_reflection_functions()

            #Calculate the gammap parameter
            self.gamma0=self._find_gamma()

            #Accelerate the calculation of the albedo
            self._accelerate_lambertian_albedo()

    def get_albedo(self,eta,zeta,delta,lamb,**params):
        return self._get_albedo(eta)
        
    #####################################
    #Complimentary routines
    #####################################
    def _load_reflection_functions(self):
        """Load value of reflection fucntions.

        Update:
            fint: 2d interpolating function:
                x: eta (cosine incident angle)
                y: zeta (cosine scattering angle)

        Notes:
            Tab. (2.3) in Sobolev (1975).
        """
        data_ss=np.loadtxt(Misc.get_data("diffuse_reflection_function.data"))
        eta=data_ss[1:,0]
        gamma=data_ss[0,1:]
        f=data_ss[1:,1:]
        self.fint=interp2d(gamma,eta,f)  

    def _calc_reflection_coefficient(self,eta,zeta,gamma0=1):
        """Reflection coefficient of a semi-infinite (tau = infinity) atmosphere with (gray) 
        single scattering albedo gamma0

        Requires:
            - _loadReflectionFunctions

        Notes:
            Ec. (2.43) in Sobolev (1975).
        """
        rho0=gamma0*self.fint(gamma0,eta)[0]*self.fint(gamma0,zeta)[0]/(4*(eta+zeta))
        return rho0

    def _calc_spherical_albedo(self,gamma0):
        """
        Compute spherical albedo from single scattering albedo for a semi-infinite atmosphere.

        Parameters:
            gamma0: single scattering albedo (0<=gamma0<=1), float.

        Returns:
            AS: ratio of the energy diffusely reflected by a spherical planet (0<=AS<=1), float.

        Requires:
            - _loadReflectionFunctions

        Notes:
            Ec. (1.87) in Sobolev (1975).    
        """

        AS=4*dblquad(lambda y,x,*args:self._calc_reflection_coefficient(x,y,*args)*x*y,
                     0,1,lambda x:0,lambda x:1,epsrel=1e-2,args=(gamma0,))[0]
        return AS

    def _find_gamma(self):
        """
        Starting with a target spherical albedo AS, find the value of the single scattering albedo gamma0
        of a semi-infinite atmosphere having that Albedo.

        Returns:
            gamma0: the value of gamma0 corresponding to AS (0<=gamma0<=1), float.
        """
        if np.isclose(self.AS,1,rtol=1e-2):
            return 1
        function=lambda gamma0:self._calc_spherical_albedo(gamma0)-self.AS
        gamma0=bisect(function,0.0,1.0,rtol=1e-4)
        return gamma0 if gamma0<=1 else 1

    def _calc_lambertian_albedo(self,eta):
        """
        Notes: 
            Yanovistkii (1973)
        """
        integrand=lambda zeta:self._calc_reflection_coefficient(eta,zeta,gamma0=self.gamma0)*zeta
        AL=2*quad(integrand,0,1,epsrel=1e-3)[0]
        return AL

    def _accelerate_lambertian_albedo(self):
        etas=np.linspace(0,1,20)
        ALs=np.array([self._calc_lambertian_albedo(eta) for eta in etas])
        self._get_albedo=interp1d(etas,ALs)

