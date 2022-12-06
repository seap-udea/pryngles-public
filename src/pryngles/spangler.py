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

#Aliases
sci=Science

import pandas as pd
import random

#Specialized plotting methods
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib import animation
from celluloid import Camera # getting the camera
import itertools
from tqdm import tqdm



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Spangler
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Spangler(PrynglesCommon):
    """A Spangler associated to an object or set of objects.
        
       There are two ways to initialize a Spangler:
        
            Creating a Spangler for a single object:
            
                Mandatory:
    
                    nspangles: int, default = 0:
                        Number of spangles in spangling.
    
                Optional:
    
                    body_hash: string, default = None:
                        Hash identifying the body to which spangles are associated 
                        (see Body documentation for explanation about hash).
    
                    spangle_type: int, default = 0:
                        Type of spangle (see *_SPANGLE in Consts module).
    
                    n_equ: numpy Array (3), default = [0,0,1]:
                        unitary vector normal to {equ} (equatorial) plane.
    
                    alpha_equ: float, default = 0:
                        Roll angle of x-axis of equatorial system (not implemented yet)
    
                    center_equ: numpy Array (3), default = [0,0,0]:
                        Position of the spnagler in the {equ} (equatorial) system.
    
                    center_ecl: numpy Array (3), default = [0,0,0]:
                        Position of the spnagler in the {ecl} (ecliptic) system.
                        
                    w, q0: float [rad/ut, rad], default = 0, 0:
                        Angular velocity and reference latitude at t = 0.
    
            Joining a set of Spanglers (several objects):
    
                spanglers: list of Spanglers. default = []:
                    Set of spanglers to join.
    
    Core attributes:
    
        nspangles: int:
            Total number of spangles.
    
        data: Pandas DataFrame: 
            Dataframe containing all the information about the spangling.
            For Columns see global variable SPANGLER_COLUMNS.
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    def __init__(self,
                 #Initialization using specific options
                 #Initialization with a list of spanglers
                     spanglers=[],
                 #Basic
                     nspangles=1,
                     name=None,
                     n_equ=SPANGLER_COLUMNS["n_equ"],
                     alpha_equ=SPANGLER_COLUMNS["alpha_equ"],
                     center_equ=SPANGLER_COLUMNS["center_equ"],
                     center_ecl=SPANGLER_COLUMNS["center_ecl"],
                 #Optional
                     w=SPANGLER_COLUMNS["w"],
                     q0=SPANGLER_COLUMNS["q0"],
                ):
        
        #Common attributes
        self.n_obs=np.array([0,0,1])
        self.n_luz=np.array([0,0,1])
        self.d_obs=self.d_luz=1
        self.sample=None
        self.geometry=-1
        
        #Direction of vantages point in spherical coordinates
        self.rqf_obs=sci.spherical(self.n_obs)
        self.rqf_luz=sci.spherical(self.n_luz)
        self.center_luz=None
        self.center_obs=None
        
        #Transformation matrices from equatorial to ecliptic coordinates
        self.M_equ2ecl=dict()
        self.M_ecl2equ=dict()
        
        #Convex hulls of spanglers
        self.qhulls=dict()
        
        #Required for plotting
        self.fig2d=None
        self.ax2d=None
        self.fig3d=None
        self.ax3d=None
        
        #Create a spanglers with a list of other spanglers
        if len(spanglers)>0:
            verbose(VERB_SIMPLE,f"Joining {len(spanglers)} spanglers")
            self._join_spanglers(spanglers)
            
        #Create a spangler with the desired options
        else:
            #Attributes
            self.nspangles=nspangles
            self.shape="vanilla" #No geometry defined
            
            #Default property values
            self._defaults=deepcopy(SPANGLER_COLUMNS)

            if not name:
                #Generate a random hash for object
                self.name=str(random.getrandbits(16))
                verbose(VERB_VERIFY,f"Generating random hash {self.name}")
            else:
                self.name=name
                
            self._defaults.update(dict(name=self.name))
            
            #Update other parameters
            self._defaults.update(
                dict(w=w,q0=q0)
            )

            #Create Spangler dataframe
            if self.nspangles>0:
                
                #Create a simple DataFrame with the default values
                self.data=pd.DataFrame([list(self._defaults.values())]*self.nspangles,
                                       columns=self._defaults.keys())

                #Update positions
                self.set_positions(
                    n_equ=n_equ,alpha_equ=alpha_equ,
                    center_equ=center_equ,center_ecl=center_ecl,
                    t=None
                )
        
            else:        
                verbose(VERB_SIMPLE,f"Creating a blank Spangler")
                #Creat a blank DataFrame
                self.data=pd.DataFrame(columns=self._defaults.keys())
        
    def reset_state(self):
        """Reset spangler state
        """
        self.data[list(SPANGLER_SOURCE_STATES)+list(SPANGLER_VISIBILITY_STATES)]=False
        self.data["unset"]=True
        for coords in "int","obs","luz":
            self.data["hidden_by_"+coords]=""
            self.data["transit_over_"+coords]=""

    def set_scale(self,scale):
        """Set scale
        
        scale: float:
            New scale.  All lengths will be multiplied by scale, areas by scale**2 and
            vector components by scale.
        """
        self.scale=scale
        self.data[SPANGLER_LENGTHS]*=self.scale
        self.data[SPANGLER_AREAS]*=self.scale**2
        for vector in SPANGLER_VECTORS:
            self.data[vector]=[np.array(v)*scale for v in self.data[vector]]
        
    def _join_spanglers(self,spanglers):
        """
        Join spanglers into a single spangler

        Parameters:
            spanglers: list of Spanglers:
                Spanglers to join.
        """
        self.name=[]
        for spangler in spanglers:
            if not isinstance(spangler,Spangler):
                raise AssertionError(f"One of the spangler is not an Spangler instance")
                
            if spangler.name in self.name:
                raise ValueError(f"Hash '{spangler.name}' already included in spangler '{self.name}'")
                
            self.name+=[spangler.name]

        #When joining there is no single geometry
        self.shape="Join"
        
        #Set of spanglers
        self.spanglers=spanglers

        #Concatenate data
        datas=[spangler.data for spangler in spanglers]
        self.data=pd.concat(datas,ignore_index=True)

        self.M_equ2ecl=dict()
        for spangler in spanglers:
            self.M_equ2ecl.update(spangler.M_equ2ecl)

        #Join properties
        self.nspangles=len(self.data)
        
    def get_mem_usage(self,info=False):
        """Get size of the Spangler DataFrame in megabytes

        Optional parameters:
            info: boolean, default = False:
                Get detailed information on memory usage.
                
        Return:
            mem_usage: float:
                Size in Mb of the memomory usage.
        """
        mem_usage=self.data.memory_usage(deep=True).sum()/1024**2
        if info:
            print(f"Basic information:\n")
            self.data.info(memory_usage="deep")
            print(f"\nDetailed size:\n")
            print(self.data.memory_usage(deep=True).to_string())
            print(f"\nTotal size: {mem_usage:.1g} Mb")

        return mem_usage
        
    def set_positions(self,
                      n_equ=[],alpha_equ=0,
                      center_equ=[],center_ecl=[],
                      t=None
                     ):
        """
        Set the positions and orientation of spanglers in all reference systems.

        Parameters:

            n_equ: list/array (3), default = []:
                Normal vector towards north pole equatorial system.

            alpha_equ: float, default = 0:
                Roll angle of x-axis of equatorial system (not implemented yet)

            center_equ: list/array (3), default = []:
                Location of the center of the body with respect to the barycenter in the equatorial system.

            center_ecl: list/array (3), default = []:
                Location of the center of the body with respect to the barycenter in the ecliptic system.

            t: float, default = None:
                Time.  This quantity is used to update the equatorial coordinates.
                If None, equatorial coordinates are not set.

        Return:
            None

        Update:
    
            If n_equ:
                Rotation matrices M_equ2ecl

            If t is provided:
                Coordinates of the spangles in the equatorial, (x_equ,y_equ,z_equ).
                Normals to the spangle (ns_equ)

            In all cases:
                Coordinates of the spangles, (x_ecl,y_ecl,z_ecl).
            
        """
        verbose(VERB_VERIFY,f"Setting positions")

        #Update normal vectors
        qupdate=False

        #Update center
        if len(center_equ)>0:
            verbose(VERB_VERIFY,f"Updating center in {{equ}} to {center_equ}")
            self.data["center_equ"]=[center_equ]*self.nspangles
            
        if len(center_ecl)>0:
            verbose(VERB_VERIFY,f"Updating center {{ecl}} to {center_ecl}")
            self.data["center_ecl"]=[center_ecl]*self.nspangles

        if len(n_equ)>0:
            verbose(VERB_VERIFY,f"Generating equatorial transformation matrices from n_equ = {n_equ}")

            #Unitary equatorial vector
            n_equ,one=spy.unorm(n_equ)
            self.data["n_equ"]=[n_equ]*self.nspangles

            #Transformation matrices
            self.M_equ2ecl[self.name],M_ecl2equ=sci.rotation_matrix(n_equ,alpha_equ)

            qupdate=True

        #Update equatorial coordinates by rotation
        if t is not None:
            verbose(VERB_VERIFY,f"Updating rotations at t = {t}")

            self.data["q_equ"]=[q+q0+w*t for q,w,q0 in zip(self.data.q_equ,self.data.w,self.data.q0)]
            self.data[["x_equ","y_equ","z_equ"]]=                [sci.cartesian(r) for r in np.array(self.data[["r_equ","q_equ","f_equ"]])]

            qupdate=True

        #If equatorial positions have been changed
        if qupdate and self.sample:
 
            #Update spangles orientations
            verbose(VERB_VERIFY,f"Generating normal vectors")

            #If the spangler has been poputaled update normals
            if self.sample:
                self.data["ns_equ"]=pd.Series(
                    list(
                        self.sample.update_normals(self.data[["x_equ","y_equ","z_equ"]])
                    ),dtype=object
                )

        #Convert from equatorial to ecliptic
        verbose(VERB_VERIFY,f"Converting to equatorial")
        self.data[["x_ecl","y_ecl","z_ecl"]]=            [np.matmul(self.M_equ2ecl[sph],r+cequ)+cecl             for sph,r,cequ,cecl in zip(self.data.name,
                                        np.array(self.data[["x_equ","y_equ","z_equ"]]),
                                        self.data.center_equ,self.data.center_ecl)]
        
        #Update orientation of the spangle
        self.data["ns_ecl"]=[np.matmul(self.M_equ2ecl[sph],n) for sph,n in zip(self.data.name,
                                                                               self.data.ns_equ)]
        
        #Update matrix of the transformation from ecliptic to local (horizontal) reference frame of the spangle

        #Search all spangles pointing towards ez or -ez
        cond=pd.Series([abs(spy.vdot(ns,[0,0,1]))!=1 for ns in self.data.ns_ecl])
        if sum(cond)>0:
            verbose(VERB_VERIFY,f"Setting local vectors based on ns: {sum(cond)}")
            #wy = ez x ns because with this definition seen from above the system is oriented as usual
            self.data.loc[cond,"wy_ecl"]=pd.Series([spy.unorm(spy.vcrss([0,0,1],ns))[0]                                                     for ns in self.data[cond].ns_ecl],dtype=object).values
            self.data.loc[cond,"wx_ecl"]=pd.Series([spy.vcrss(wy,ns)                                                     for ns,wy in zip(self.data[cond].ns_ecl,self.data[cond].wy_ecl)],
                                                   dtype=object).values
            cond=~cond
        else:
            cond=[True]*self.nspangles

        #Spangles pointing towards ez or -ez
        if sum(cond)>0:
            verbose(VERB_VERIFY,f"Setting local matrix based on ex: {sum(cond)}")
            self.data.loc[cond,"wx_ecl"]=pd.Series([[1,0,0]]*sum(cond),dtype=object).values
            self.data.loc[cond,"wy_ecl"]=pd.Series([spy.unorm(spy.vcrss(ns,wx))[0]                                                     for ns,wx in zip(self.data[cond].ns_ecl,self.data[cond].wx_ecl)],
                                                   dtype=object).values
            
        #Update velocities
        #Not implemented yet

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Tested methods from module file spangler
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def populate_spangler(self,
                          shape="circle",preset=False,spangle_type=SPANGLE_SOLID_ROCK,
                          scale=1,seed=0,**shape_args):
        
        """Populate data of a Spangler using points generated with a given geometry.
        
        Parameters:
                
            shape: string, default = "circle":
                Shape of the Sampler.  Available: "circle", "ring", "sphere".
    
            spangle_type: int, default = SPANGLE_SOLID_ROCK:
                Type of spangle.  See Constants module for a list of spangle types.
    
            preset: boolean, default = False:
                If true the spangler is populated with preset data (see class Sampler for details).
                
            scale: float. default = 1:
                Scale size of the object.
                
            seed: integer. default = 0:
                Value of the integer seed of random number generation (if 0 no random seed is set).
                If a non-zero seed is used the position of the spangle will be always the same.
                
            shape_args: dictionary:
                See Sampler methods documentation.
                 
        """
        #Check if preset
        if preset:
            verbose(VERB_VERIFY,f"Populating spangler from preset for {shape}")
            preset=(shape,shape_args)
            self.sample=Sampler(preset=preset,N=self.nspangles,seed=seed)   
        else:
            verbose(VERB_VERIFY,f"Generating spangler from scratch")
            self.sample=Sampler(N=self.nspangles,seed=seed)
            exec(f"self.sample.gen_{shape}(**shape_args)")
    
        self.shape=shape
        self.data["geometry"]=self.sample.geometry
        self.data["spangle_type"]=spangle_type
    
        if self.sample.geometry in [SAMPLER_GEOMETRY_SPHERE]:
            #Purge sample if it is in 3d
            verbose(VERB_VERIFY,f"Purging 3d sample")
            self.sample.purge_sample()
            self.nhidden=0
            
        elif self.shape == "ring":
            #Number of hidden points
            """
            The number of hidden points for a ring is choosen in such a way that the ratio between t
            he area of the circle sector to the area of the circle segment is larger than (1-epsilon)
            
            Ag / As = (r^2 sin(teta)/2)/(r^2 teta/2) > (1-epsilon)
            
            Expanding sin(teta) we get:
            
            1 - teta^2/6 > 1 - epsilon
            """
            nhidden=int(2*np.pi/(6*SPANGLER_EPS_BORDER)**0.5)
    
            #Add hidden spangles to ring inner borders
            pp_border=np.zeros((nhidden,3))
            ss_border=np.zeros((nhidden,3))
            ns_border=np.zeros((nhidden,3))
            for i,theta in enumerate(np.linspace(0,2*np.pi,nhidden)):
                pp_border[i]=[self.sample.ri,theta,0]
                ss_border[i]=[self.sample.ri*np.cos(theta),
                              self.sample.ri*np.sin(theta),
                              0]
                ns_border[i]=[0,0,1]
            self.sample.pp=np.vstack((self.sample.pp,pp_border))
            self.sample.ss=np.vstack((self.sample.ss,ss_border))
            self.sample.ns=np.vstack((self.sample.ns,ns_border))
            self.sample.N+=nhidden
            self.nhidden=nhidden
                    
        #Check if number of samples is not equal to that of spangles defined when the spangler was created
        if self.sample.N!=self.nspangles:
            verbose(VERB_SYSTEM,f"Sample size {self.sample.N} is different from spangles {self.nspangles}. Adjusting.")
    
            #Difference between sampler number of points and number of spanglers
            dif=self.sample.N-self.nspangles
            
            if dif>0:
                #Add spangles
                verbose(VERB_SYSTEM,f"Adding {dif} entries to DataFrame")
                for i in range(dif):
                    df=pd.DataFrame([self.data.iloc[-1]])
                    self.data=pd.concat([self.data,df],ignore_index=True)
            else:
                #Remove spangles
                verbose(VERB_SYSTEM,f"Removing {-dif} entries to DataFrame")
                self.data.drop(range(self.nspangles+dif,self.nspangles),inplace=True)
                
            self.nspangles=self.sample.N
        
        #Area
        self.data["asp"]=self.sample.aes*scale**2
        self.data["dsp"]=2*(self.data["asp"]/np.pi)**0.5
        
        #Update scale
        self.data["scale"]=scale
    
        #Store positions in DataFrame
        self.data[["x_equ","y_equ","z_equ"]]=self.sample.ss*scale
        self.data[["r_equ","q_equ","f_equ"]]=self.sample.pp
        self.data["r_equ"]*=scale
    
        #Update normal vectors
        self.data["ns_equ"]=pd.Series(list(self.sample.ns),dtype=object)
            
        #Hide border points in case of ring
        if shape == "ring":
            self.data.loc[self.nspangles-self.nhidden:self.nspangles,"hidden"]=True
            
        #Update positions
        self.set_positions()
        
    def _update_column_order(self):
        """Reorder columns in a more convenient way.
        """
        for col in self.data.columns:
            if col not in SPANGLER_KEY_ORDERING:
                raise AssertionError(f"Column {col} not present in key ordering.  Please checl SPANGLER_KEY_ORDERING")
        
        #Order columns
        self.data=self.data.loc[:,SPANGLER_KEY_ORDERING]
         
    
    def plot3d(self,
               coords="ecl",
               only=None,
               center_at=None,
               not_plot=[],
               fsize=5,
               factor=1.2,
               statemark=0,
               show_directions=False
              ):
        """Plot spangle in 3d.
    
        Optional parameters:
        
            coords: list of strings, default = ["x_ecl","y_ecl","z_ecl"]:
                Which coordinates do you want to plot.  
                Available: equ, ecl, obs, luz, int.
                
            only: string, default = None:
                Plot only the object with this hash.        
    
            center_at: string, default = None:
                Hash of the object around which the plotting will be centered at (see name column
                of the Spangler DataFrame).
                
            not_plot: list of strings, default = []:
                List of object hashes to not plot.
                
            fsize: tuple (2), default = 5:
                Size of the figure.  The parameter figsize used at creating the figure will be 
                figsize = (fsize,fsize).
    
            factor: float, default = 1.2:
                Size of the coordinate axes.  factor = 1 correspond to axis equal to maximum and minumum.
                
            statemark: float, default = 0:
                If different than 0 mark with state the spangles in 3d plot.  
                It will mark the 1-markstate spangles in plot.
                
            show_directions: boolean, default = False:
                If True show the direction of normal vectors to spangles and vectors directed from the 
                origin of the intersection system of reference to the spangles.
                
        Color coding:
            
            Determinative of color:
            
                By default or in darkness: color of darkness (dark blue)
            
                If illuminated: color of the spangle.
            
                If in shadow: color of shadow.
            
            Modification of the color: 
            
                If not visible: reduce level of color to half
    
        """
        bgcolor='k'
    
        #Plot only a given object
        if only:
            not_plot=list(self.data.name.unique())
            if only not in not_plot:
                raise ValueError(f"Spangler '{only}' not among available spanglers ({not_plot})")
            else:
                not_plot.remove(only)
                center_at=only
    
        #Check if plot is in the ecliptic system
        qecl=True
        if 'ecl' not in coords:
            qecl=False
    
        scoords=coords
        coords=[f"x_{scoords}",f"y_{scoords}",f"z_{scoords}"]
        
        #Center
        cond=(self.data.name==center_at)
        x_cen,y_cen,z_cen=self.data[cond][coords].mean() if sum(cond)>0 else np.array([0,0,0])
        
        #Figure
        fig=plt.figure(figsize=(fsize,fsize))
        fig.patch.set_facecolor(bgcolor)
        ax=fig.add_subplot(111,projection='3d',facecolor=bgcolor)
        ax.axis("off")
        
        #Spangles
        for i in range(self.nspangles):
    
            #Avoid plotting 
            name=self.data.loc[i,"name"]
            if name in not_plot:
                continue
            
            #Reference transparency of spangles
            alpha_base=0.5
    
            #Avoid hidden spangles
            if self.data.loc[i,"hidden"]:
                continue
    
            spangle_type=self.data.loc[i,"spangle_type"]
    
            #Define the color according to illumination or shadow
            state=""
            color_hls=SPANGLES_DARKNESS_COLOR #Default color: gray
    
            #Define color according to illumination or shadow
            if self.data.loc[i,"shadow"]:
                #Inside a shadow
                state+="S."
                color_hls=SHADOW_COLOR_LUZ #Gray
            
            if self.data.loc[i,"illuminated"]:
                #Illuminated
                state+="I."
                color_hls=SPANGLE_COLORS[spangle_type] #Planet color
                
            #Modify color according to visibility, transmission or darkness
            if not self.data.loc[i,"illuminated"]:
                #In darkness
                state+="D."
                color=Plot.rgb(color_hls) #No color modification
                
            if self.data.loc[i,"transmit"]:
                #Transmitting
                state+="T."
                color=Plot.rgb(color_hls) #No color modification
    
            if not self.data.loc[i,"visible"]:
                #Not visible
                state+="N."
                color=Plot.rgb([color_hls[0],
                                color_hls[1]/2, #Reduce level to half
                                color_hls[2]
                               ])
            else:
                #Invisible
                state+="V."
                color=Plot.rgb(color_hls) #No color modification
                
            if self.data.loc[i,"unset"]:
                state+="U."
                color=Plot.rgb([0,0.5,0])
                
            #Define alpha according to albedo
            alpha=alpha_base*self.data.albedo_gray_normal[i]
    
            center=[self.data[coords[0]][i]-x_cen,self.data[coords[1]][i]-y_cen,self.data[coords[2]][i]-z_cen]
            radius=self.data.dsp[i]/2
            zDir=self.data[f"ns_{scoords}"][i]
    
            #verbose(VERB_DEEP,i,center,radius,zDir)
            Plot.circle3d(ax,
                          center=center,
                          radius=radius,
                          zDir=zDir,
                          color=color,alpha=alpha,lw=0)
            if statemark:
                if np.random.rand()>1-statemark:
                    ax.text(center[0],center[1],center[2],state,fontsize=6,color='w')
            
        #Aspect
        ax.set_box_aspect([1,1,1])
    
        #Zoom around center
        cond=(self.data.name==center_at)
        cond=cond if sum(cond)>0 else [True]*self.nspangles
    
        #Not 
        cond=cond&(~self.data.name.isin(not_plot))
        
        #Range
        maxval=1.0*np.abs(self.data[cond][coords].to_numpy()-[x_cen,y_cen,z_cen]).max()
        ax.set_xlim(-maxval,maxval)
        ax.set_ylim(-maxval,maxval)
        ax.set_zlim(-maxval,maxval)
        
        #Decoration
        xmin,xmax=factor*np.array(list(ax.get_xlim()))
        ymin,ymax=factor*np.array(list(ax.get_ylim()))
        zmin,zmax=factor*np.array(list(ax.get_zlim()))
    
        #Axis
        ax.plot([xmin,xmax],[0,0],[0,0],'w-',alpha=0.3)
        ax.plot([0,0],[ymin,ymax],[0,0],'w-',alpha=0.3)
        ax.plot([0,0],[0,0],[zmin,zmax],'w-',alpha=0.3)
        ax.text(xmax,0,0,rf"$x_{{{scoords}}}$",color='w',alpha=0.5,fontsize=8)
        ax.text(0,ymax,0,rf"$y_{{{scoords}}}$",color='w',alpha=0.5,fontsize=8)
        ax.text(0,0,zmax,rf"$z_{{{scoords}}}$",color='w',alpha=0.5,fontsize=8)
        
        #Plot n_obs and n_luz vector only in the case of ecliptic system
        increase=1.05*factor*maxval
        if qecl:
            ax.quiver(+self.n_luz[0]*increase,+self.n_luz[1]*increase,+self.n_luz[2]*increase,
                      -self.n_luz[0]*increase,-self.n_luz[1]*increase,-self.n_luz[2]*increase,
                      color='y',alpha=0.7)
            ax.text(self.n_luz[0]*increase,self.n_luz[1]*increase,self.n_luz[2]*increase,
                    r"$n_{luz}$",color='w',alpha=0.7,fontsize=8,ha='left',va='bottom')
            ax.quiver(+self.n_obs[0]*increase,+self.n_obs[1]*increase,+self.n_obs[2]*increase,
                      -self.n_obs[0]*increase,-self.n_obs[1]*increase,-self.n_obs[2]*increase,
                      color='c',alpha=0.7)        
            ax.text(self.n_obs[0]*increase,self.n_obs[1]*increase,self.n_obs[2]*increase,
                    r"$n_{obs}$",color='c',alpha=0.7,fontsize=8,ha='right',va='top')
            ax.view_init(30,60)
        else:
            r_obs,t_obs,f_obs=sci.spherical(self.n_obs)
            ax.view_init(f_obs*Consts.rad,t_obs*Consts.rad)
            
        #Show vectors
        if show_directions:
            
            cond=cond&(~self.data.hidden)
            
            if scoords=="ecl":
                nstr="n_int_"+scoords
            else:
                nstr="n_"+scoords
            
            #It is important to stress that "vectors" should be directed towards the light-source 
            vectors=np.array(list(self.data.loc[cond,nstr]))
            
            normals=np.array(list(self.data.loc[cond,"ns_"+scoords]))
            
            ax.scatter(self.data.loc[cond,"x_"+scoords],self.data.loc[cond,"y_"+scoords],self.data.loc[cond,"z_"+scoords],
                       marker="*",c='w',s=10)
            ax.quiver(self.data.loc[cond,"x_"+scoords],self.data.loc[cond,"y_"+scoords],self.data.loc[cond,"z_"+scoords],
                      vectors[:,0],vectors[:,1],vectors[:,2],color='w',label="From intersection")
            ax.quiver(self.data.loc[cond,"x_"+scoords],self.data.loc[cond,"y_"+scoords],self.data.loc[cond,"z_"+scoords],
                      normals[:,0],normals[:,1],normals[:,2],color='y',label="Normal")
            
            leg=ax.legend(loc='lower right',facecolor='k',ncol=3,prop={'size':8},bbox_to_anchor=(0.5,-0.01, 0.5, 0.5))
            frame=leg.get_frame()
            frame.set_edgecolor("k")
            for text in leg.get_texts():
                text.set_color("w")
    
        #Title
        ax.set_title(f"Spangler {self.shape}, N = {self.nspangles}",
                     color='w',fontsize=10)
        Plot.pryngles_mark(ax)
        
        fmark=""
        if statemark:
            fmark=f", I/D: Illum./Dark, V/N: Visible/Invisible, S: Shadow, T: Transmit"
        
        #Scale
        ax.text2D(0,0,f"Axis scale: {maxval*factor:.2g} {fmark}",
                fontsize=7,color='w',
                transform=ax.transAxes)
    
        fig.tight_layout()
        self.fig3d=fig
        self.ax3d=ax
    
    def set_intersect(self,
                      nvec=[0,0,1],
                      alpha=0,
                      center=None,
                      name=None,
                     ):
        """Set the positions and orientation of spanglers in an intersection direction
    
        Parameters:
    
            nvec: list/array (3), default = [0,0,1]:
                Vector pointing towards the vantage point from where the intersection will be computed. 
                It can be normalized or not.  The components are in the ecliptic reference system.
                            
            alpha: float, default = 0:
                Roll angle of x-axis.
                
            center: list/array (3), default = None:
                Location of the vantage point in the ecliptic reference system.
                If None, we assume that the vantage point is at an infinite distance.
                
            name: string, default = None:
                Spangler hash to which the transformation will be applied.
    
        Return:
        
            cond: boolean array:
                Over which spangles the transformation was applied.
                
            n_int: array (3):
                Normal vector towards the vantage point.
                
            d_int: float:
                Distance to vantage point.  If 'center' is None, this distance is set to numpy.inf.
                
        Create:
        
            qhulls: dictionary:
                Convex hulls of bodies from this vantage point.
                
                key: 
                    name
                
                value: 
                    list with hulls corresponding to each name.
                
        Update:
    
            Coordinates of the spangles in the intersection system, (x_int,y_int,z_int).
    
            Normal to spangles in the intersection system, ns_int.
            
            
        Notes:
            If the intersection direction is in the center of the body (for instance, when a ring or a bubble 
            is illuminated from the center), set intersect to True for all spangles and compute the distance and
            relative orientation (cos_int) of the spangles correspondingly.
    
        """
        
        verbose(VERB_SIMPLE,
                f"Setting intersect using nvec = {nvec}, alpha = {alpha} center = {center}, name = {name}")
        
        verbose(VERB_VERIFY,f"Generating intersection matrices from pvec = {nvec}")
    
        #Unitary observer vector
        n_int,norm=spy.unorm(nvec)
        alpha_int=alpha
        
        #Store n_int and d_int for update state purposes
        self.rqf_int=sci.spherical(n_int)
        self.n_int=n_int
        
        #Distance to center of intersection
        if center is None:
            self.infinite=True
            d_int=np.inf
            center=np.array([0,0,0])
        else:
            self.infinite=False
            d_int=np.linalg.norm(center)
            center=np.array(center)
        self.d_int=d_int
    
        #Transformation matrices
        M_int2ecl,self.M_ecl2int=Science.rotation_matrix(n_int,alpha_int)
        self.M_int2ecl=M_int2ecl
        
        #Depending on body
        cond=[True]*self.nspangles
        if name:
            cond=(self.data.name==name)
    
        #If no point is of type name
        if sum(cond)==0:
            return
            
        #Update positions in the intersection reference frame
        self.data.loc[cond,["x_int","y_int","z_int"]]=        [np.matmul(self.M_ecl2int,r-center) for r in np.array(self.data[cond][["x_ecl","y_ecl","z_ecl"]])]
        
        #Center of the object in the observer reference system
        center_int=[np.matmul(self.M_ecl2int,c_ecl+np.matmul(self.M_equ2ecl[sp],c_equ)-center)               for sp,c_ecl,c_equ in zip(self.data[cond].name,
                                             np.array(self.data[cond].center_ecl),
                                             np.array(self.data[cond].center_equ))]
        self.data.loc[cond,"center_int"]=pd.Series(center_int).values
        
        #According to distance to intersetcion point generate z_cen_int
        if self.infinite:
            self.data.loc[cond,"z_cen_int"]=-np.inf
        else:
            self.data.loc[cond,"z_cen_int"]=np.array(center_int)[:,2]
    
        #Pseudo-cylindrical coordinates in the observer system
        self.data.loc[cond,["rho_int","az_int","cosf_int"]]=        [sci.pcylindrical(r) for r in          np.array(self.data[cond][["x_int","y_int","z_int"]])-np.vstack(self.data[cond].center_int)]
    
        #Compute distance to intersection of each spangle and the 
        if self.infinite:
            #Distance to all points is assumed infinite
            self.data.loc[cond,"n_int"]=pd.Series([[0,0,1]]*sum(cond),dtype=object)
            self.data.loc[cond,"d_int"]=np.inf
            self.data.loc[cond,"n_int_ecl"]=pd.Series([list(n_int)]*sum(cond),dtype=object)
        else:
            #Distance to origin of coordinates in the int system where the center is located
            lista=[spy.unorm(list(-r)) for r in np.array(self.data[cond][["x_int","y_int","z_int"]])]
            self.data.loc[cond,["n_int","d_int"]]=pd.DataFrame(lista,columns=["n_int","d_int"],
                                                               index=self.data[cond].index)
            n_int_ecl=[spy.mxv(M_int2ecl,n_int) for n_int in self.data.n_int[cond]]
            self.data.loc[cond,"n_int_ecl"]=pd.Series(n_int_ecl,dtype=object)
            
        #Azimuth of the direction of the intersection vector in the tangent plane of the spangle
        self.data.loc[cond,"azim_int"]=[np.arctan2(spy.vdot(wy,n_int),
                                                   spy.vdot(wx,n_int)) \
                                        for wy,wx,n_int in zip(self.data.wy_ecl[cond],
                                                               self.data.wx_ecl[cond],
                                                               self.data.n_int_ecl[cond])]
        
        #Update spangles orientations
        lista=[np.matmul(self.M_ecl2int,n_ecl) for n_ecl in self.data[cond].ns_ecl]
        self.data.loc[cond,"ns_int"]=pd.Series(lista,dtype=object).values
        
        #Cosine of the direction of the intersection vector and the normal to the spangle
        if self.infinite:
            #In this case n_int is a global variable
            self.data.loc[cond,"cos_int"]=[spy.vdot(n_ecl,n_int) for n_ecl in self.data.ns_ecl[cond]]
        else:
            #In this case n_int is a per-spangle variable
            self.data.loc[cond,"cos_int"]=[np.vdot(ns,n_int)                                        for ns,n_int in zip(self.data.ns_int[cond],self.data.n_int[cond])]
            
        #Set areas
        self.data.loc[cond,"asp_int"]=self.data.loc[cond,"asp"]
        
        return cond,n_int,d_int
    
    def _calc_qhulls(self):
        
        """Compute convex hulls for a given intersection configuration
        """
        
        #Convex hulls
        for name in Misc.flatten([self.name]):
    
            self.qhulls[name]=[]
            cond_obj=(self.data.name==name)
            center=list(self.data[cond_obj].center_int.iloc[0])
            zord=min(self.data[cond_obj].z_int)
    
            if (self.data[cond_obj].hidden).sum()==0:
    
                #Convex hull of whole objects
                cond_hull=(cond_obj)&(~self.data[cond_obj].hidden)
                verbose(VERB_SIMPLE,"Hull points (whole object):",sum(cond_hull))
                qhull=Science.get_convexhull(self.data[cond_hull][["x_int","y_int"]])
                vhull=qhull.volume if qhull else 0
                
                self.qhulls[name]+=[dict(
                    name=name,
                    hulltype="cen",
                    center=center,
                    zord=zord,
                    qhull=qhull,
                    vhull=vhull,
                )]
    
            else:
                #Convex hull of objects with a hole (eg. rings)
    
                #Plane of rings
                cond_hidden=(cond_obj)&(self.data[cond_obj].hidden)
                hidden=self.data[cond_hidden][["x_int","y_int","z_int"]].values
                nhidden=len(hidden)
                p1,p2,p3=hidden[0],hidden[int(nhidden/3)],hidden[2*int(nhidden/3)]
                plane=Plane(p1,p2,p3)
    
                #Convex hull of hidden points (the hole)
                cond_hull=(cond_obj)&(self.data[cond_obj].hidden)
                verbose(VERB_SIMPLE,"Hull points (hidden):",sum(cond_hull))
                qhull=Science.get_convexhull(self.data[cond_hull][["x_int","y_int"]])
                vhull=qhull.volume if qhull else 0
    
                self.qhulls[name]+=[dict(
                    name=name,
                    hulltype="hidden",
                    center=center,
                    zord=zord,
                    qhull=qhull,
                    vhull=vhull,
                    plane=plane
                )]
    
                #Convex hull of no hidden points
                cond_hull=(cond_obj)&(~self.data[cond_obj].hidden)
                verbose(VERB_SIMPLE,"Hull points (visible ring):",sum(cond_hull))
                qhull=Science.get_convexhull(self.data[cond_hull][["x_int","y_int"]])
                vhull=qhull.volume if qhull else 0
    
                self.qhulls[name]+=[dict(
                    name=name,
                    hulltype="plane",
                    center=center,
                    zord=zord,
                    qhull=qhull,
                    vhull=vhull,
                    plane=plane
                )]
               
    
    def set_observer(self,nvec=[0,0,1],alpha=0,center=None):
        """Set the positions and orientation of spanglers in the observer system.
    
        Parameters:
    
            nvec: list/array (3), default = [0,0,1]:
                Normal vector towards the observer.
    
            alpha: float, default = 0:
                Roll angle of x-axis of observer system.
                
            center: list/array(3), default = None:
                Define the position of the vantage point in the ecliptic system.
                
        """
        verbose(VERB_SIMPLE,f"Setting observer")
        
        #Set observer
        cond,self.n_obs,self.d_obs=self.set_intersect(nvec,alpha,center)
        
        #Set properties
        self.alpha_obs=alpha
        self.rqf_obs=sci.spherical(self.n_obs)
        self.center_obs=center.copy() if center else center
        
        self.data.loc[cond,"visible"]=False
        #self.data.loc[cond,SPANGLER_COL_OBS]=self.data.loc[cond,SPANGLER_COL_INT].values
        #"""
        self.data.loc[cond,SPANGLER_COL_OBS]=pd.DataFrame(self.data.loc[cond,SPANGLER_COL_INT].values,
                                                          columns=SPANGLER_COL_OBS,
                                                          index=self.data[cond].index)
        #"""
        #self.data.loc[cond,SPANGLER_COL_OBS]=self.data.loc[cond,SPANGLER_COL_INT]
        
        
        #Update states
        self.data.unset=False
        
        #Condition for visibility
        """
        & ! Hidden
        & z_cen_obs+scale < 0: spangle is observable from the observer vantage point-
            (
                | cos_obs > 0: spangle it is towards the observer
                | Spangle type is semitransparent
            )
        """
        cond=    (~self.data.hidden)&    ((self.data.z_cen_obs+self.data.scale)<0)&    (        (self.data.cos_obs>0)|        (self.data.spangle_type.isin(SPANGLES_SEMITRANSPARENT))
        )
        self.data.loc[cond,"visible"]=True
        
    def set_luz(self,nvec=[0,0,1],alpha=0,center=None,name=None):
        """Set the positions and orientation of spanglers in the light-source system.
    
        Parameters:
    
            nvec: list/array (3), default = [0,0,1]:
                Normal vector towards the observer.
    
            alpha: float, default = 0:
                Roll angle of x-axis of observer system.
                
            center: list/array(3), default = None:
                Define the position of the vantage point in the ecliptic system.
    
            name: string, default = None:
                Body to apply this light direction
                
        Update:
            This method update the 'illuminated' and 'transmit' states.
                    
        Note:
            For updating the 'transmit' state it is required that the observer be set.
            
        """
       
        verbose(VERB_SIMPLE,f"Setting light-source")
     
        #Set intersect of all points in order to prepare the update luz
        cond,self.n_luz,self.d_luz=self.set_intersect(nvec,alpha,center,name=None) 
        verbose(VERB_SIMPLE,f"Number of points: {sum(cond)}")
        
        #Depending on body choose which spangles to change
        cond=[True]*self.nspangles
        if name:
            cond=(self.data.name==name)
        
        #Set the light source direction in spherical coordinates
        self.rqf_luz=sci.spherical(self.n_luz)
        
        #Set the default value of the states to change in False
        self.data.loc[cond,"illuminated"]=False
        self.data.loc[cond,"transmit"]=False
        
        #Conditions
        #self.data.loc[cond,SPANGLER_COL_LUZ]=deepcopy(self.data.loc[cond,SPANGLER_COL_INT].values)
        self.data.loc[cond,SPANGLER_COL_LUZ]=pd.DataFrame(self.data.loc[cond,SPANGLER_COL_INT].values,
                                                          columns=SPANGLER_COL_LUZ,
                                                          index=self.data[cond].index)
        
        #Set relative azimuth
        self.data.loc[cond,"azim_obs_luz"]=self.data.loc[cond,"azim_obs"]-self.data.loc[cond,"azim_luz"]
        
        #Update states
        self.data.loc[cond,"unset"]=False
        
        #Condition for illumination
        """
        & ! Hidden
        & z_cen_luz+scale < 0: spangle is in front of the light-source.
            (
                | geometry = circle : 2d spangles are always illuminated
                | spangle_type = stellar: stellar spangles are always illuminated
                | cos_luz > 0: spangle it is towards the light source
            )
        """
        cond=    cond&    (~self.data.hidden)&    ((self.data.z_cen_luz+self.data.scale)<0)&    (        (self.data.geometry==SAMPLER_GEOMETRY_CIRCLE)|        (self.data.cos_luz>0)|        (self.data.spangle_type==SPANGLE_STELLAR)|        (self.data.spangle_type.isin(SPANGLES_SEMITRANSPARENT))
        )
        self.data.loc[cond,"illuminated"]=True
    
        #Conditions for transmission:
        """
        & No hidden
        (
            & Spangle type is semitransparent
            & cos_obs . cos_luz < 0: observer and light source are in opposite sides
        )
        
        ATTENTION: TRANSMISSION IS ONLY PROPERLY SET IF OBSERVER HAVE BEEN PREVIOUSLY SET.
        """
        cond=    cond&    (~self.data.hidden)&    (     (self.data.spangle_type.isin(SPANGLES_SEMITRANSPARENT))&     ((self.data.cos_luz*self.data.cos_obs)<=0)
        )
        self.data.loc[cond,"transmit"]=True
        
    
    def plot2d(self,
               coords="obs",
               center_at=None,
               include=[],
               exclude=[],
               axis=True,
               fsize=5,
               newfig=True,
               show_azim=False,
               highlight=None,
               maxval=None
              ):
        """
        Plot spangle.
    
        Basic parameters:
        
            coords: string, default = obs:
                which coordinates do you want to use. Available: 'equ', 'ecl', 'int', 'obs', 'luz'.
        
        Other parameters:
                
            center_at: string, default = None:
                Hash of the object around which the plotting will be centered at (see name column
                of the Spangler DataFrame).
                
            include: string, default = None:
                List of objects (hashes) to plot exclusively.
                        
            exclude: list of strings, default = []:
                List of objects (hashes) to not plot.
                
            fsize: integer, default = 5:
                Size of figure
                
            newfig: boolean, default = True:
                If True a new figure is created.  The False value is intended for animations.
                
            show_azim: boolean, default = False:
                If True show azimuth of the observer and light source direction on each spangle.
                
            highlight: tuple, default = None:
            
                A tuple containing:
        
                    1) A boolean mask telling which spangles to highlight
                    2) A dictionary with options for a scatter command
                    
                Example:
                    
                    #Highlight all the spangles belonging to star which are visible
                    cond=(sys.data.name=="Star")&(sys.data.visible)
                    sys.sg.plot2d(highlight=(cond,dict(c='c')))
                    
                If the dictionary is empty it uses the default dict(s=1,c='w')
                
            maxval: float, default = None:
                Change range of plot.  If None it is calculated automatically.
        """
        
        #Global properties of the plot
        bgcolor='k'
        fig_factor=fsize/5
    
        #Create figure and axes
        if "fig2d" not in self.__dict__ or newfig:
            fig=plt.figure(figsize=(fsize,fsize))
            fig.patch.set_facecolor(bgcolor)
            ax=fig.add_subplot(111,facecolor=bgcolor)
    
            #Keep figure and axe
            self.fig2d=fig
            self.ax2d=ax
    
        #Convert list 
        include_string=[]
        for body in include:
            if isinstance(body,PrynglesCommon):
                include_string+=[body.name]
            else:
                include_string+=[body]
        include=include_string
    
        exclude_string=[]
        for body in exclude:
            if isinstance(body,PrynglesCommon):
                exclude_string+=[body.name]
            else:
                exclude_string+=[body]
        exclude=exclude_string
        
        if isinstance(center_at,PrynglesCommon):
            center_at=center_at.name
    
        #Plot only a set of objects
        if len(include)>0:
            
            #List of spanglers names
            exclude=list(self.data.name.unique())
            
            for name in include:
                if name not in exclude:
                    raise ValueError(f"Spangler '{name}' not among available spanglers ({exclude})")
                else:
                    exclude.remove(name)
                    
            #Center at the first object in the list
            center_at=include[0]
    
        #Center of plot
        cond=(self.data.name==center_at)
        x_cen,y_cen,z_cen=self.data[cond][[f"x_{coords}",f"y_{coords}",f"z_{coords}"]].mean() if sum(cond)>0 else np.array([0,0,0])
    
        #Select plotting bodies
        cond_included=(~self.data.hidden)&(~self.data.name.isin(exclude))
        num_included=sum(cond_included)
        if num_included==0:
            raise AssertionError(f"No body remain after removing {exclude}")
        data=self.data[cond_included]
        
        #Calculate range of plot
        cond_maxval=(~data.hidden)&(~data.name.isin(exclude))
        cond_maxval=cond_maxval if sum(cond_maxval)>0 else [True]*num_included        
        if not maxval:
            maxval=1.2*np.abs(np.array(data[cond_maxval][[f"x_{coords}",f"y_{coords}"]])-np.array([x_cen,y_cen])).max()
    
        #Function to determine the size of the spangles
        size_factor=1/2.5
        size_points=lambda dsp,cos_obs:size_factor*(dsp[cond])*abs(cos_obs)**0.5
        
        ##########################################################
        #Plotting properties according to state
        ##########################################################
        #Default colors and sizes
        colors=np.array(['#000000']*num_included)
        sizes=np.array([0.0]*num_included)
        marker='o'
        #All
        """
        cond=[True]*len(data)
        sizes[cond]=0
        """
    
        #Illuminated
        cond=(data.visible)&(data.illuminated)
        verbose(VERB_SIMPLE,f"Visible and illuminated: {cond.sum()}")
        colors[cond]=[Plot.rgb([SPANGLE_COLORS[stype][0],
                                SPANGLE_COLORS[stype][1]*min((cos_luz*cos_obs+0.3),1),
                                SPANGLE_COLORS[stype][2]],
                                to_hex=True) for stype,cos_luz,cos_obs in zip(data[cond].spangle_type,
                                                                           abs(data[cond].cos_luz),
                                                                           abs(data[cond].cos_obs))
                     ] #Object color
        sizes[cond]=size_points(data.dsp[cond],data.cos_obs[cond])
    
        #Not illuminated
        cond=(data.visible)&(~data.illuminated)
        verbose(VERB_SIMPLE,f"Visible and not illuminated: {cond.sum()}")
        colors[cond]=Plot.rgb(SPANGLES_DARKNESS_COLOR,to_hex=True)
        sizes[cond]=size_points(data.dsp[cond],data.cos_obs[cond])
        
        #In shadow
        cond=(data.visible)&(data.shadow)
        verbose(VERB_SIMPLE,f"Visible and not illuminated: {cond.sum()}")
        colors[cond]=Plot.rgb(SHADOW_COLOR_LUZ,to_hex=True)
        sizes[cond]=size_points(data.dsp[cond],data.cos_obs[cond])
    
        if coords!="obs":
            #Not visible
            cond=(~data.visible)&(data[f"z_{coords}"]>0)
            colors[cond]=Plot.rgb(SHADOW_COLOR_OBS,to_hex=True)
            sizes[cond]=size_points(data.dsp[cond],data.cos_obs[cond])
    
        #Transmitting
        cond=(data.visible)&(data.transmit)&(data.illuminated)
        verbose(VERB_SIMPLE,f"Visible, illuminated and transmitting: {cond.sum()}")
        colors[cond]=[Plot.rgb([SPANGLE_COLORS[stype][0],
                                SPANGLE_COLORS[stype][1]*min((cos_luz*cos_obs+0.3),1)/2,
                                SPANGLE_COLORS[stype][2]],
                                to_hex=True) for stype,cos_luz,cos_obs in zip(data[cond].spangle_type,
                                                                           abs(data[cond].cos_luz),
                                                                           abs(data[cond].cos_obs))
                     ] #Object color
        sizes[cond]=size_points(data.dsp[cond],data.cos_obs[cond])
        
        ##########################################################
        #Dot sizes scaled according to figure size
        ##########################################################
        ppd=72./self.ax2d.figure.dpi
        trans=self.ax2d.transData.transform
        dot_size=lambda x:int(((trans((1,x/maxval))-trans((0,0)))*ppd)[1]**2)
        
        ##########################################################
        #Plotting properties according to state
        ##########################################################
        st=[max(dot_size(s),0.1) if s>0 else 0 for s in sizes]
        sargs=dict(c=colors,s=st,marker=marker,zorder=-100)
        self.ax2d.scatter(data[f"x_{coords}"]-x_cen,data[f"y_{coords}"]-y_cen,**sargs)
        
        #Ranges
        self.ax2d.set_xlim(-maxval,maxval)
        self.ax2d.set_ylim(-maxval,maxval)
        
        factor=1
        xmin,xmax=factor*np.array(list(self.ax2d.get_xlim()))
        ymin,ymax=factor*np.array(list(self.ax2d.get_ylim()))
        
        ##########################################################
        #Show azim
        ##########################################################
        if show_azim:
            #Choose which directions to show
            cond=(self.data["cos_"+coords]>=0)&(~self.data.hidden)&(cond_included)&(self.data["spangle_type"]!=SPANGLE_STELLAR)
    
            #Options of arrows showing direction
            quiver_args=dict(scale=15,scale_units='width',
                             width=0.005,alpha=0.6,zorder=+1000,headwidth=0)
            
            #Quiver plot of azimuth for light
            azx=[mh.cos(x) for x in self.data[cond].azim_luz]
            azy=[mh.sin(x) for x in self.data[cond].azim_luz]
            
            
            self.ax2d.quiver(self.data[cond]["x_"+coords]-x_cen,self.data[cond]["y_"+coords]-y_cen,
                             azx,azy,color='m',label="Az.luz",**quiver_args)
    
            #Quiver plot of azimuth for observer
            azx=[mh.cos(x) for x in self.data[cond].azim_luz]
            azy=[mh.sin(x) for x in self.data[cond].azim_luz]
            self.ax2d.quiver(self.data[cond]["x_"+coords]-x_cen,self.data[cond]["y_"+coords]-y_cen,
                             azx,azy,color='w',label="Az.obs",**quiver_args)
    
            #Quiver plot of elevation for light
            tx=np.sqrt(1-self.data[cond].cos_int**2).values
            ty=self.data[cond].cos_int.values
            self.ax2d.quiver(self.data[cond]["x_"+coords]-x_cen,self.data[cond]["y_"+coords]-y_cen,
                             tx,ty,color='c',label="Elev.luz",**quiver_args)
    
            #Legend decoration
            leg=self.ax2d.legend(loc='lower right',facecolor='k',ncol=3,prop={'size':8},
                                 bbox_to_anchor=(0.5, -0.05, 0.5, 0.5))
            frame=leg.get_frame()
            frame.set_edgecolor("k")
            for text in leg.get_texts():
                text.set_color("w")
            axis=False
    
        ##########################################################
        #Highlight spangles
        ##########################################################
        if highlight:
            if len(highlight)<2:
                raise AssertionError("Highlight should include conditions and scatter options")
            
            def_args_scatter=dict(c='w',s=0.1,marker='*')
            cond_highlight,args_scatter=highlight
            def_args_scatter.update(args_scatter)
            self.ax2d.scatter(self.data[cond_highlight&cond_included]["x_"+coords]-x_cen,
                              self.data[cond_highlight&cond_included]["y_"+coords]-y_cen,
                              **def_args_scatter)
            
        ##########################################################
        #Show axis and other labels
        ##########################################################
        if newfig and axis:
            self.ax2d.plot([xmin,xmax],[0,0],'w-',alpha=0.3)
            self.ax2d.plot([0,0],[ymin,ymax],'w-',alpha=0.3)
            self.ax2d.text(xmax,0,fr"$x_{{{coords}}}$",color='w',alpha=0.5,fontsize=8*fig_factor)
            self.ax2d.text(0,ymax,fr"$y_{{{coords}}}$",color='w',alpha=0.5,fontsize=8*fig_factor)
    
            #Scale
            center_text=""
            if center_at:
                center_text=f", Center at '{center_at}'"
            self.ax2d.text(0,0,f"Axis scale: {maxval*factor:.2g}{center_text}",
                      fontsize=8*fig_factor,color='w',
                      transform=self.ax2d.transAxes)
    
        ##########################################################
        #Decorate plot
        ##########################################################
        if newfig:
            #Title
            label_obs=""
            lamb=0
            phi=0
            if coords=="obs":
                lamb=self.rqf_obs[1]*Consts.rad
                phi=self.rqf_obs[2]*Consts.rad
                coords_label=f"($\lambda$,$\\beta$) : ({lamb:.1f}$^\circ$,{phi:.1f}$^\circ$)"
            elif coords=="luz":
                lamb=self.rqf_luz[1]*Consts.rad
                phi=self.rqf_luz[2]*Consts.rad
                coords_label=f"($\lambda$,$\\beta$) : ({lamb:.1f}$^\circ$,{phi:.1f}$^\circ$)"
            elif coords=="int":
                lamb=self.rqf_int[1]*Consts.rad
                phi=self.rqf_int[2]*Consts.rad
            coords_label=f"($\lambda$,$\\beta$) : ({lamb:.1f}$^\circ$,{phi:.1f}$^\circ$)"
    
            if coords=="ecl":
                coords_label=""
    
            label_obs=f"{coords} {coords_label}"
            self.ax2d.text(0.5,1.01,f"{label_obs}",
                         transform=self.ax2d.transAxes,ha='center',
                         color='w',fontsize=10*fig_factor)
    
            self.ax2d.axis("off")
            Plot.pryngles_mark(self.ax2d)
        
        ##########################################################
        #Adjust sizes
        ##########################################################
        self.ax2d.axis("equal")
        self.fig2d.tight_layout()
        
        return x_cen,y_cen
    
    
    def update_intersection_state(self,excluded=[],included=[]):
        """Update state of intersections
        
        Otional Parameters:
            exluded: list, default = []:
                List of objects to exclude from the calculation.
                
            included: list, default = []:
                Objects to include in the calculation.
        
        """    
        #Update qhulls using the latest intersection state
        self._calc_qhulls()
        
        #Check if at least one qhull has been computed
        if len(self.qhulls) == 0:
            raise AssertionError("You must set an intersection vantage point.")
    
        #List of objects in spangler
        names=list(Misc.flatten([self.name]))
        
        if len(included):
            excluded=[n for n in names if n not in included]
        verbose(VERB_SIMPLE,f"Exclusion list: {excluded}")
        
        #Objects included when performing intersection calculation
        cond_included=self.data.name.apply(lambda x:x not in excluded)
        if cond_included.sum()==0:
            raise AssertionError("You have excluded all objects when calculating intersetions.")
        else:
            verbose(VERB_SIMPLE,f"Points included in calculation: {cond_included.sum()}")
        
        #Under the current circumstances all this spangles are intersecting 
        cond=    (~self.data.hidden)&    (     (self.data.cos_int>0)|     (self.data.spangle_type.isin(SPANGLES_SEMITRANSPARENT))    )&    (cond_included)
        
        self.data.loc[cond,"intersect"]=True
        self.data.hidden_by_int=""
        self.data.transit_over_int=""
            
        #Loop over objects producing intersection
        for name in names:
            
            #If body is excluded
            if name in excluded:
                verbose(VERB_SIMPLE,f"Skipping {name} from intersection computation")
                continue
            
            #Points in present body
            cond=(self.data.name==name)
            geometry=self.data[cond].geometry.iloc[0]
            scale=self.data[cond].scale.iloc[0]
            
            #If this body is not in the field-of-view, avoid computation
            z_cen_int=self.data[cond].z_cen_int.iloc[0]
            scale=self.data[cond].scale.iloc[0]
            if (z_cen_int+scale)>=0:
                continue
    
            #By default for all objects
            inhull_not_in_hole=[True]
            
            verbose(VERB_SIMPLE,f"Calculating intersections for '{name}'")
            
            #Hull area of the body
            ahull=0
            for i,hull in enumerate(self.qhulls[name]):
                
                qhull=hull["qhull"]
                if qhull is None:
                    verbose(VERB_SIMPLE,f"No hull for '{name}'")
                    continue
                vhull=hull["vhull"]
                
                htype=hull["hulltype"]
                xcen,ycen,zcen=hull["center"]
                zord=hull["zord"]
                
                verbose(VERB_SIMPLE,f"Hull {i+1} for '{name}' of type '{htype}'")
    
                #Evaluate conditions
                inhull=sci.points_in_hull(self.data[["x_int","y_int"]],qhull)&(~cond)&(cond_included)
                below=np.array([False]*self.nspangles)
                above=np.array([False]*self.nspangles)
                
                if htype=="hidden":
    
                    #Holes
                    inhull_not_in_hole=(~inhull)
                    verbose(VERB_SIMPLE,f"Points outside hidden hull for '{name}': {sum(inhull_not_in_hole)}")
                    hull["notinhole"]=sum(inhull_not_in_hole)
    
                    continue
                    
                else:
                    #Body
                    verbose(VERB_SIMPLE,f"Points not in hole for '{name}:{htype}': {sum(inhull_not_in_hole)}")
        
                    #Area of the body
                    ahull+=vhull
                    
                    #Spangles to evaluate
                    cond_vis=(self.data.cos_int>0)|(self.data.spangle_type.isin(SPANGLES_SEMITRANSPARENT))
                    cond_int=(~self.data.hidden)&(self.data.name!=name)&(cond_vis)
    
                    if htype=="cen":
                        below=(inhull_not_in_hole)&(inhull)&(self.data[cond_int]["z_int"]<=zcen)
                        above=(inhull_not_in_hole)&(inhull)&(self.data[cond_int]["z_int"]>zcen)
                        
                    elif htype=="plane":
                        #Not in hole, inhull, not hidden, not in object and intersect
                        cond_full=(inhull_not_in_hole)&(inhull)&(cond_int)
                        verbose(VERB_SIMPLE,"Fulfilling all conditions:",sum(cond_full))
                        
                        plane=hull["plane"]
                        below[cond_full]=[plane.is_below(r,[0,0,1]) for r in self.data[cond_full][["x_int","y_int","z_int"]].values]
                        above[cond_full]=~below[cond_full]
                    else:
                        raise ValueError("Type of hull '{htype}' not recognized")
                
                #Store information
                verbose(VERB_SIMPLE,f"Points in hull for '{name}:{htype}': {sum(inhull)}")
                verbose(VERB_SIMPLE,f"Points below '{name}:{htype}': {sum(below)}")
                
                #Set visibility
                self.data.loc[below,"intersect"]=False
                self.data.loc[below,"hidden_by_int"]=self.data.loc[below,"hidden_by_int"]+f"{name}:{zord:.3e}&"
                    
                #Compute distance to center for transiting spangles
                self.data.loc[above,"string_int"]=[f"{name}:{zord:.3e}:{((r[0]-xcen)**2+(r[1]-ycen)**2)**0.5/scale:.3e}&"                                                for r in self.data[above][["x_int","y_int"]].values]
                self.data.loc[above,"transit_over_int"]=self.data.loc[above,"transit_over_int"]+self.data.loc[above,"string_int"]
                
                hull["inhull"]=sum(inhull)
                hull["below"]=sum(below)
                
            #Correct areas
            if geometry != SAMPLER_GEOMETRY_CIRCLE:
                cond=(self.data.name==name)&(self.data.cos_int>=0)&(~self.data.hidden)
                
                #This what the sum actually is
                ahull_expected=(self.data.loc[cond,"asp_int"]*abs(self.data.loc[cond,"cos_int"])).sum()
                #This is the expected value
                ahull=np.pi*scale**2
                
                #Normalization factor to get sum = ahull
                norma=ahull/ahull_expected
                
                #Final area
                self.data.loc[cond,"asp_int"]*=norma
        
    def update_visibility_state(self):
        """Update states and variables related to visibility
        """
        self.update_intersection_state()
        self.data[SPANGLER_COL_OBS]=self.data[SPANGLER_COL_INT]
        
        self.data.hidden_by_obs=self.data.hidden_by_obs+self.data.hidden_by_int
        self.data.transit_over_obs=self.data.transit_over_obs+self.data.transit_over_int
        self.data.visible=self.data.visible&self.data.intersect
    
    def update_illumination_state(self,excluded=[],included=[]):
        """Update states and variables related to illumination
        """
        self.update_intersection_state(excluded,included)
        self.data[SPANGLER_COL_LUZ]=self.data[SPANGLER_COL_INT]
        
        self.data.hidden_by_luz=self.data.hidden_by_luz+self.data.hidden_by_int
        self.data.transit_over_luz=self.data.transit_over_luz+self.data.transit_over_int
        
        #Update illumination
        self.data.illuminated=self.data.illuminated&self.data.intersect
        
        #Not intersected and spangles in the direction of the light-source are for sure shadowed spangles
        #Condition for visibility
        """
        | shadow
        & ~intersect : spangle does not intersect
            (
                | cos_obs > 0: spangle it is towards the observer
                | Spangle type is semitransparent
            )
        """
        self.data.shadow=    (self.data.shadow)|    (
            (~self.data.intersect)&\
            (
                (self.data.cos_int>0)|
                (self.data.spangle_type.isin(SPANGLES_SEMITRANSPARENT))
            )
        )
        
        #Stellar spangles are always illuminated
        cond=(self.data.spangle_type==SPANGLE_STELLAR)
        self.data.loc[cond,"illuminated"]=True
        
        #In stellar spangles cos_luz = cos_obs for not having strange visual representations
        self.data.loc[cond,"cos_luz"]=self.data.loc[cond,"cos_obs"]
    
