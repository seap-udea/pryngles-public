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

import rebound as rb
from tqdm import tqdm


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class System
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class System(PrynglesCommon):
    """Creates a planetary system.
    
        Initialization attributes:
    
            units: list of strings, default = ['au','msun','yr2pi']:
                Units used in calculations following the conventions and signs of rebound.
                The order SHOULD always be MKS: length, mass, time (in that order)
    
        Optional attributes:
    
            resetable: boolean, default = False:
                If True the system is resetable, namely you can reset it to the initial system.
                
            filename: string, default = None:
                File to load system.
    
        Derived attributes:
    
            sim: Class Simulation:
                Rebound Simulation object.
    
            ul, um, ut: float [SI units]:
                Value of the conversion factors for each unit.
    
            G: float [ul^3/ut^2/um]
                Value of the gravitational constant.
    
            bodies: dictionary:
                Bodies in the system.
    
            nbodies: int:
                Number of bodies.
    
            nparticles: int:
                Numbre of particles in rebound simulation.
    
            spangler: Class Spangler:
                Spangler object with all the spangles in the system.
    
        Examples:
    
            #Create a system
            sys=System(units=["au","msun","yr"])
            sys.sim.integrator='whfast'
            sys.sim.dt=0.01
    
            #Add star (by default, m = 1)
            S=sys.add()
    
            #Add planet, when an object is added, it is automatically spangled
            P=sys.add("Planet",radius=0.1,m=1e-3,a=1,e=0.2)
    
            #Add moon: orbital elements are respect to ecliptic system
            M=sys.add("Planet",parent=P,radius=0.01,m=1e-7,a=0.1,e=0.01)
    
            #Add ring system
            R=sys.add("Ring",parent=P,fi=1.5,fe=2.5,albedo_gray_normal=0.5,tau_gray_optical=3)
    
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def __init__(self,
                 filename=None,
                 units=['au','msun','yr2pi'],
                 resetable=False
                ):
        
        if filename:
            self.load_from(filename)
            return
        
        #Rebound simulation
        self.sim=None
        self._simulated=False
        
        #Attributes by default
        
        #List of bodies in the system
        self.bodies=odict()
        
        #Root of the tree of bodies
        self.root=None
        
        #Center of the light-source in the system
        self.source=None
        self.center_root=np.array([0,0,0])
        
        #Orbital configuration
        self.orbital_configuration=None
        
        #Observer properties
        self.n_obs=np.array([0,0,1])
        self.alpha_obs=0  
        self.center_obs=None
        
        #Check if spangled
        self._spangled=False
        
        #Check if observer has been set
        self._observer_set=False
        self._luz_set=False
        
        #Initialize spangler object
        self.sg=None
        
        #Is the system resetable?
        self._resetable=resetable
        if self._resetable:
            #Create temporary file
            self._snap_file_name = "/tmp/pryngles-system.pkl"
        
        #Update rebound units
        self.update_units(units)
        
        #By default spangle scatterers
        """
        This is the list of the class of scatterers used to calculate the scattering in
        different types of spangles.
        
        The structure of this dictionary is:
        
            key: integer (or enumerator):
                This is the column spangle_type in the spangler DataFrame.
                
            item: tuple (2):
                Component 1: 
                    Class of scatterer.
                Component 2: 
                    Dictionary mapping the initialization properties of the scatterer to columns
                    in the spangler DataFrame.
        
        Example of item:
        
            SPANGLE_ATMOSPHERIC:(LambertianGrayAtmosphere,dict(AS="albedo_gray_spherical"))
            
                This means that for spangles of the type SPANGLE_ATMOSPHERIC Pryngles will 
                instantiate an object of the class LambertianGrayAtmosphere.  This class have a
                single parameter, the spherical albedo AS.  The dictionary means that when 
                instantiating the object the column "albedo_gray_spherical" will be used to 
                initialize the object.                
        """
        self.spangle_scatterers={
            SPANGLE_ATMOSPHERIC:(LambertianGrayAtmosphere,dict(AS="albedo_gray_spherical")),
            SPANGLE_GRANULAR:(LambertianGraySurface,dict(AL="albedo_gray_normal")),
            SPANGLE_LIQUID:(LambertianGraySurface,dict(AL="albedo_gray_normal")),
            SPANGLE_SOLID_ICE:(LambertianGraySurface,dict(AL="albedo_gray_normal")),
            SPANGLE_SOLID_ROCK:(LambertianGraySurface,dict(AL="albedo_gray_normal")),
            SPANGLE_GASEOUS:(BlackBodySurface,dict()),
            SPANGLE_STELLAR:(BlackBodySurface,dict()),
        }
        
    def update_units(self,units):
        """Update units of the system
        """
        #Check units
        if units[0] not in rb.units.lengths_SI:
            raise ValueError(f"Length unit provided '{units[0]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.lengths_SI.keys())}")
        if units[1] not in rb.units.masses_SI:
            raise ValueError(f"Mass unit provided '{units[1]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.masses_SI.keys())}")
        if units[2] not in rb.units.times_SI:
            raise ValueError(f"Time unit provided '{units[2]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.times_SI.keys())}")
        
        #Units        
        self.units=units
        self._ul,self._um,self._ut=self.units
        #self.sim.units=self.units
        
        #Canonical units of the system
        self.ul=rb.units.convert_length(1,self._ul,"m")
        self.um=rb.units.convert_mass(1,self._um,"kg")

        #Compute the units of time
        sim=rb.Simulation()
        sim.units=self.units
        self.G=sim.G
        self.ut=np.sqrt(self.G*self.ul**3/(self.um*GSI))
        
        #Update system
        self._update_system()
        
    def _get_source(self,body):
        """Get the source of light (stellar body) in the center of a body
        """
        if (body.parent is None) or (body.kind == "Star"):
            return body

        elif body.parent.kind == "Star":
            return body.parent

        else:
            return self._get_source(body.parent)

    def _update_system(self):
        """Update system properties
        """
        self.nbodies=len(self.bodies)
        if self._simulated:
            self.nparticles=len(self.sim.particles)
        
    def _is_spangled(self):
        """Check if system is spangled
        """
        return True if self.sg else False
    
    def reset_state(self):
        """Reset the state of the spangler
        """
        self.sg.reset_state()
        self._observer_set=False
        self._luz_set=False

    def save_to(self,filename):
        """Save system from file
        
        Parameters:
            filename: string:
                Path to file where the object will be pickled.
                
        Result:
            File 'filename' for regular object and 'filename.rbin' for rebound simulation
        """
        if self._simulated:
            #Rebound file
            rb_filename=filename+".rbin"

            #Save rebound state
            verbose(VERB_SIMPLE,"Saving rebound simulation")
            self.sim.save(rb_filename)

            #Since rebound have ctypes it cannot be pickled
            del self.sim
            self._simulated=True

        #Pickle system
        PrynglesCommon.save_to(self,filename)

        if self._simulated:
            #Load again rebound
            self.sim=rb.Simulation(rb_filename)

    def load_from(self,filename):
        """Load system from filename
                
        Parameters:
            filename: string:
                Path to file where the object will be pickled.
                There to be 2 files: 'filename' (with the regular object) and filename.rbin with 
                rebound simulation.
        """
        #Load system
        self=PrynglesCommon.load_from(self,filename)

        if self._simulated:
            #Rebound file
            rb_filename=filename+".rbin"

            #Load rebound
            verbose(VERB_SIMPLE,"Loading rebound simulation")
            self.sim=rb.Simulation(rb_filename)
        
    def status(self):
        if self._simulated:
            print(f"System with {self.nbodies} bodies and {self.nparticles} particles (rings and disk are not particles)")
            self.sim.status()
        else:
            print(f"Simulation for this system has not been yet initialized. Use System.initialize_simulation()")


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Tested methods from module file scatterer
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def update_scatterers(self):
        """Update the scatterers of the spangles
        """
        if not self._spangled:
            raise AssertionError("You need to spangle the system before updating the scatterers.")
        
        #Update scatterer only for the non-assigned one
        cond=(self.data.scatterer=="")
        for index in self.data[cond].index:
            #Get spangle
            spangle=self.data.loc[index]
            #Get spangle sype
            spangle_type=spangle["spangle_type"]
            #Get scatterer class and options description
            spangle_scatterer,spangle_options=self.spangle_scatterers[spangle_type]
            #Build options of scatterers from options description
            scatterer_options={**dict(zip(spangle_options.keys(),spangle[list(spangle_options.values())]))}
            #Instantiate object of scatterer and save hash into DataFrame
            self.data.loc[index,"scatterer"]=spangle_scatterer(**scatterer_options).hash
    

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Tested methods from module file system
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def add(self,kind="Star",parent=None,**props):
        """Add an object to the system
        
        Examples:
        
            sys=System()
            S=sys.add("Star",m=2)
        
        Parameters:
        
            kind: string, default = "Star":
                Kind of object: Star, Planet, Ring (see BODY_KINDS).
        
            parent: Body, default = None:
                Parent object of the body.
                
            props: dictionary:
                List of properties of the body.
                
        Returns:
            
            Body
                Body added to the system.
                
        Examples:
            #Add star (by default, m = 1)
            S=sys.add()
    
            #Add planet, when an object is added, it is automatically spangled
            P=sys.add("Planet",radius=0.1,m=1e-3,x=1,vy=0.2)
    
            #Add moon: orbital elements are respect to ecliptic system
            M=sys.add("Planet",parent=P,radius=0.01,m=1e-7,a=0.1,e=0.01)
    
            #Add ring system
            R=sys.add("Ring",parent=P,fi=1.5,fe=2.5,albedo_gray_normal=0.5,tau_gray_optical=3)        
            
        """
        if kind is None:
            raise AssertionError("You must provide a valid object kind (Star, Planet, Ring).")
    
        if kind not in BODY_KINDS:
            raise ValueError(f"Object kind '{kind}' is not recognized.")
    
        #Legacy
        if 'primary' in props:
            parent=props["primary"]
        if kind=="Observer":
            parent=self.root
        
        #Default parameters
        if self.root:
            if (kind!="Star") and (parent is None):
                parent=self.root
                
            if kind=="Planet":
                if "m" not in props:
                    props["m"]=0.1*parent.m
                if "radius" not in props:
                    props["radius"]=0.5*parent.radius
                if "a" not in props:
                    if "a" in parent.__dict__:
                        props["a"]=0.5*parent.a
                    
        #Create body
        props.update(dict(name_by_kind=True))
        self.__body=eval(f"{kind}(parent=parent,**props)")
        
        if self.__body.name in self.bodies:
            raise ValueError(f"An object with name '{self.__body.name}' has been already added.")
        
        #If we have a root object and no parent has been provided
        if not parent:
            if self.root:
                raise ValueError(f"A root object alread exist in the system ({self.root.name}) and you do not provided a parent body for {self.__body.name}.")
            else:
                self.root=self.__body
                verbose(VERB_SIMPLE,f"Setting the root object as {self.root.name}")
            
        self.bodies[self.__body.name]=self.__body
        
        #Create the shined body tree
        self.__body.shined=[]
    
        #Update system
        self._update_system()
        
        #Set the source of the object
        if self.__body.source:
            #Check that the source is a body
            if not isinstance(self.__body.source,Body):
                raise ValueError(f"The source of body must be an actual Body.")        
            #Check that the source is among the bodies
            if self.__body.source.name not in self.bodies:
                raise ValueError(f"The source of body {self.__body.name} is not among system bodies {list(self.bodies.keys())}.")
            #Check that the source be a star
            if self.__body.source.kind!="Star":
                raise ValueError(f"The source of body {self.__body.name} must be a Star.  You set {self.__body.source.name} which is a {self.__body.source.kind}.")
            self.__body.source.shined+=[self.__body.name]
        else:
            if self.__body.kind=="Star":
                self.__body.source=self.__body
            elif self.__body.parent==self.root:
                self.__body.source=self.root
            else:
                self.__body.source=self.__body.parent.source
        self.__body.source.shined+=[self.__body.name]
        
        verbose(VERB_SIMPLE,f"Object '{kind}' with name '{self.__body.name}' has been added.")
        return self.__body
    
    
    def initialize_simulation(self,orbital_tree=None,**rebound_options):
        """Initialize rebound simulation using a given orbital tree.
        
        Parameters:
            orbital_tree: list of pairs, default = None:
                A scheme showing how the bodies in the system are organized as in a 
                hierarchical N-body system (see OrbitUtil.build_system).
                
                Examples:
                    Simple system: star (S), planet (P):
                        orbital_tree = [S,P]
                    
                    System with two planets: star (S), planet 1 (P1), planet 2 (P2):
                        orbital_tree = [[S,P1],P2]
                        
                    System with moon: star (S), planet (P), moon (M):
                        orbital_tree = [S,[P,M]]
                        
                    System with two planets and moons: star (S), planet 1 (P1), moon planet 1 (M), planet 2 (P2):
                        orbital_tree = [[S,[P1,M]],P2]
                        
                    System with two stars and one planet per star:
                        orbital_tree = [[S1,PS1],[S1,PS2]]
        
        Return:
            orbit: object Orbit:
                Object containing the hierarchical N-body system.
        
        Update:
            self.sim: Rebound Simulation:
                Simulation of the system.
                
        Note:
        
            You can 
        """
        
        #Compile orbital configuration
        if orbital_tree is None:
            i=0
            for name,body in odict(reversed(list(self.bodies.items()))).items():
                if body == self.root:
                    continue
                if body.kind == "Ring":
                    continue
                if i == 0:
                    self.orbital_tree=body
                else:
                    self.orbital_tree=[body,self.orbital_tree]
                i+=1
            self.orbital_tree=[self.root,self.orbital_tree]
        else:
            self.orbital_tree=orbital_tree
        
        #Set the rebound hash of all bodies
        for name,body in self.bodies.items():
            if body.kind=="Ring":
                body.rbhash=body.parent.name
            else:
                body.rbhash=body.name
        
        #Check that all bodies in system is in the orbital tree
        bodies=list(Misc.flatten(self.orbital_tree))
        for name,body in self.bodies.items():
            if body.kind=="Ring":
                continue
            if body not in bodies:
                raise AssertionError(f"Body '{name}' is in System but not in orbital tree.")
            
        #Build hierarchical N-body system
        orbit,pelements=OrbitUtil.build_system(self.orbital_tree,self.units)
        orbit.calculate_orbit()
        
        #Initialize simulation
        self.sim=rb.Simulation()
        self.sim.units=self.units
        
        #Add particles to simulation
        orbit.sim.move_to_com()
        for i,p in enumerate(orbit.sim.particles):
            self.sim.add(
                hash=bodies[i].name,
                m=bodies[i].m,
                x=p.x,y=p.y,z=p.z,
                vx=p.vx,vy=p.vy,vz=p.vz
            )
        self.sim.orbit=orbit
        self._simulated=True
        self._update_system()
        
        return orbit
    
    def remove(self,name):
        """Remove a body from a system.
    
        Parameters:
            name: string
                Hash of the body to remove
        
        Notes: 
            Remove eliminate body and all the childs and the childs of the childs.
    
        Example:
            sys=System()
            S=sys.add(m=2)
            sys.remove(name=S.name)
        """
        
        if name in self.bodies:
            verbose(VERB_SIMPLE,f"Removing object {name} from system")
    
            obj=self.bodies[name]
    
            #Get the list of child hashes before removing (it changes during for)
            child_hashes=list(obj.childs.keys())
            
            #Remove child objects
            for child_hash in child_hashes:
                if child_hash in self.bodies:
                    self.remove(child_hash)
                    
            #Remove object from Rebound simulation
            if obj.kind != "Ring":
                if self._simulated:
                    if self.nparticles:
                        verbose(VERB_SIMPLE,f"Removing particle {name} from simulation")
                        self.sim.remove(hash=name)
            
            #Remove object from childs of its parent
            if obj.parent:
                del obj.parent.childs[name]
            
            #Remove object from bodies
            del self.bodies[name]
    
            #Update system
            self._update_system()
        else:
            raise ValueError(f"No object with hash '{name}' in the system")
    
    def spangle_system(self):
        """Generate the spangles of the objects in the system
        
        Attributes created:
            
            spanglers: dictionary of Spangler objects:
                Spangler corresponding to each object in the system.
                
            sp: Spangler:
                Spangler corresponding to all system.
                
        Result:
            
            This method create the spangler of the system
    
        """
        if not self._simulated:
            raise AssertionError("Before spangling the system you must initialize the simulation: System.initialize_simulation().")
        
        self._spanglers=dict()
        
        #Add spangles
        for name,body in self.bodies.items():
            
            verbose(VERB_SIMPLE,f"Spangling body '{name}' (kind '{body.kind}')")
            body.spangle_body()
    
            #Center object around its position according to rebound
            body.center_ecl=np.array(self.sim.particles[body.rbhash].xyz)
            body.sg.set_positions(center_ecl=body.center_ecl)
            self._spanglers[name]=body.sg
            
        #Set the center of the source of light for each body
        for name,body in self.bodies.items():
            body.center_source=body.source.center_ecl
            if body==self.root:
                self.center_root=body.source.center_ecl
                
        #Join spanglers
        self.sg=Spangler(spanglers=list(self._spanglers.values()))
    
        #An usefule alias
        self.data=self.sg.data
        
        #Set default observer
        self.update_perspective(n_obs=self.n_obs,alpha_obs=self.alpha_obs)
        
        #Save state of the system
        if self._resetable:
            self.save_to(self._snap_file_name)
        
        #Already spangled
        self._spangled=True
    
    def _set_observer(self,nvec=[0,0,1],alpha=0,center=None):
        """Set the position of the observer
        """
        #Only set observer if it is spangled
        if self._is_spangled():
            
            #At changing the observer, reset state
            self.sg.reset_state()
            
            #Set observer
            self.sg.set_observer(nvec=nvec,alpha=alpha,center=center)
            
            #Update areas of the spangles
            
            #Update system properties
            self.d_obs=self.sg.d_obs
            self.n_obs=self.sg.n_obs.copy()
            self.rqf_obs=self.sg.rqf_obs.copy()
            self.alpha_obs=self.sg.alpha_obs
            self.center_obs=self.sg.center_obs
        
            #Update visibility
            self.sg.update_visibility_state()
            
            #Check that observer has been set
            self._observer_set=True
            
        else:
            raise AssertionError("You must first spangle system before setting observer direction.")
            
    def _set_luz_recursive(self,name,nluz,verbosity=VERB_SIMPLE):
        """Set light source for body and 
        """
        body=self.bodies[name]
        verbose(verbosity,f"Illuminating body {name}, with nluz = {nluz} and center = {body.center_source}")
        self.sg.set_luz(nvec=nluz,center=body.center_source,name=name)
        if body.childs:
            verbose(verbosity,f"\tObject {name} has childs!")
            for child_name in body.childs:
                verbose(verbosity,f"\tCalling recursively set_luz for {child_name}")
                self._set_luz_recursive(child_name,nluz)
        else:
            verbose(verbosity,f"\tObject {name} has no childs!")
                
    def _set_luz(self,verbosity=VERB_SIMPLE):
        """Set illumination in the system.
        
        Update:
            States: illuminated, shadow, hidden_by_luz
        """
        if self._is_spangled():
            
            if not self._observer_set:
                raise AssertionError("You must first set observer before setting light.")
            
            self.bodies_illuminated=[]
            for name,body in self.bodies.items():
              
                if body.kind == "Star":
                    verbose(verbosity,f"Body {body.name} is a star... skipping")
                    continue
                    
                if body.parent.name in self.bodies_illuminated:
                    verbose(verbosity,f"Parent body of {name}, {body.parent.name}, has been already illuminated")
                    continue
                            
                #Get center of body
                center=body.center_ecl
                        
                #Get source and center
                verbose(verbosity,f"Calculating illumination for '{name}' coming from '{body.source.name}' @ {body.center_source}")            
                nluz=body.center_source-center
                        
                if body.kind == "Ring" and body.parent.kind == "Star":
                    verbose(verbosity,f"Parent body of ring, {body.parent.name} is a star. All spangles will be illuminated")
                    self.sg.set_luz(nvec=nluz,center=body.center_source,name=name)
                    cond=(self.sg.data.name==name)
                    self.sg.data.loc[cond,"unset"]=False
                    self.sg.data.loc[cond,"illuminated"]=True
                    self.sg.data.loc[cond,"shadow"]=False
                else:                
                    verbose(verbosity,f"Illuminating body {name} and all its childs")
                    self._set_luz_recursive(name,nluz,verbosity)
                    self.sg.update_illumination_state(included=body.source.shined)
                    self.bodies_illuminated+=[name]
                    
            self._luz_set=True
        else:
            raise AssertionError("You must first spangle system before setting light.")
            
    def update_perspective(self,n_obs=None,alpha_obs=0,center_obs=None):
        """Update perspective (observer)
        """
        if n_obs is not None:
            #Update observing conditions
            self.n_obs,one=spy.unorm(n_obs)
            self.alpha_obs=alpha_obs
            self.center_obs=center_obs
    
        #Set observer
        self._set_observer(nvec=self.n_obs,alpha=self.alpha_obs,center=center_obs)
        self._set_luz()
        
    
    def update_body(self,body,**props):
        """Update properties of a body in the system
        
        Parameters:
            body: string or Body:
                Body to update
            
            props: dict:
                Dictionary with properties of the object
        """
        #Update spangling?
        if self._is_spangled():
            raise AssertionError("After spangling you cannot update the properties of the bodies.  Please rebuild the system")
    
        #Update body properties
        if isinstance(body,Body):
            body.update_body(**props)
        elif body in self.bodies:
            body=self.bodies[body]
            lkind=body.kind.lower()
            exec(f"body.update_{lkind}()")
        else:
            raise AssertionError("You are trying to update a body ({body}) which is not in the system")
            
        #Check if among props there is any property related to position
        if any(k in props for k in REBOUND_ORBITAL_PROPERTIES):
            raise ValueError(f"You cannot update an orbital property {props} without compromising the full simulation. Rebuild the system from scratch.")
    
    def reset(self):
        """Reset system to spangling state
        """
        if self._resetable:
            self.load_from(self._snap_file_name)
            pass
        else:
            print("System is not resetable. Use resetable = True when defining the System or when you spangle it.")
    
    def integrate(self,*args,**kwargs):
        """Integrate system
    
        Parameters:
            *args, **kwargs:
                Mandatory (non-keyword) arguments and optional (keyword) arguments for rebound.integrate.
            
        Update:
            Integrate using integrate rebound method.
            
            Update center of each body and set positions of the spangles.
        """
        #Time of integration
        t=args[0]
        verbose(VERB_SIMPLE,"Integrating up to {t}")
        
        if self._spangled:
            
            #Integrate
            self.sim.integrate(*args,**kwargs)
            self.sim.move_to_com()
        
            #Update positions
            for name,body in self.bodies.items():
                
                #Position of the body according
                body.center_ecl=np.array(self.sim.particles[body.rbhash].xyz)
    
                verbose(VERB_VERIFY,f"Updating center of body {name} @ {body.center_ecl}")
                cond=self.sg.data.name==name
                self.sg.data.loc[cond,"center_ecl"]=pd.Series([list(body.center_ecl)]*sum(cond),dtype=object).values
    
            #Update positions
            self.sg.set_positions()
            
        else:
            raise AssertionError("You must first spangle system before setting positions.")
    
    def ensamble_system(self,lamb=0,beta=0,**physics):
        """Ensamble Ringed Planet
        
        This class is for legacy purposes.
        """
        #Check if observer was provided
        if "Observer" in self.bodies:
            lamb=self.bodies["Observer"].lamb
            beta=self.bodies["Observer"].beta
            
        physics_defaults=deepcopy(LEGACY_PHYSICAL_PROPERTIES)
        physics_defaults.update(dict(limb_cs=self.bodies["Star"].limb_coeffs))
        physics_defaults.update(physics)
    
        #--CONSISTENCY--
        self._ringedplanet=dict(
            
            #Behavior
            behavior=dict(shadows=True),
            
            #Units
            CU=CanonicalUnits(UL=self.ul,UM=self.um),
    
            #Basic
            Rstar=self.bodies["Star"].radius,
            Rplanet=self.bodies["Planet"].radius,
    
            Rint=self.bodies["Ring"].fi,
            Rext=self.bodies["Ring"].fe,
            i=self.bodies["Ring"].i,
    
            a=self.bodies["Planet"].a,e=self.bodies["Planet"].e,
    
            #Orbit 
            Mstar=1,x=0,lambq=0,t0=0,kepler=False,
    
            #Observer
            eobs_ecl=np.array([lamb,beta]),
    
            #Sampling
            Np=self.bodies["Planet"].nspangles,
            Nr=self.bodies["Ring"].nspangles,
    
            Nb=0,Ns=30,
    
            #Physical properties
            physics=physics_defaults,
        )
        self.RP=RingedPlanet(**self._ringedplanet)
        return self.RP
    
