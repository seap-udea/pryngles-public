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

# # Pryngles module: Science

from pryngles import *

# ## External modules

import numpy as np
import math as mh
import spiceypy as spy

verbose=Verbose.print

# ## The Science class
# 
# The Science class is a class with routines intended to perform a wide diversity of mathematical, physical and astronomical calculations.

class Science(PrynglesCommon):pass

# ### Template method

def template(foo=1):
    """
    Method

    Parameters:

        foo: type [units], default = 1:
            Description.

    Return:

        fee: type [units]:
            Desctiption.

    """
    return foo

Science.template=template


# ### Cartesian to spherical

def spherical(xyz):
    """
    Transform cartesian coordinates into spherical coordinates

    Parameters:

        xyz: array (3):
            Cartesian coordinates

    Return:

        rqf: array (3):
            Spherical coordinates (r, theta, phi) where theta is azimutal angle and phi is 
            elevation (complement of polar angle).                

            Notice that this convention is different than that of regular vectorial calculus
            where spherical coordinates are (r,theta,phi), but theta is the polar angle and phi 
            the ezimutal one.

    """
    r,theta,phi=spy.reclat(np.array(xyz))
    theta=2*mh.pi+theta if theta<0 else theta

    return np.array([r,theta,phi])

def cartesian(rqf):
    """
    Transform cartesian coordinates into spherical coordinates

    Parameters:

        xyz: array (3):
            Cartesian coordinates

    Return:

        rqf: array (3):
            Spherical coordinates (r, theta, phi) where theta is azimutal angle and phi is 
            elevation (complement of polar angle).                

            Notice that this convention is different than that of regular vectorial calculus
            where spherical coordinates are (r,theta,phi), but theta is the polar angle and phi 
            the ezimutal one.

    """
    return spy.latrec(rqf[0],rqf[1],rqf[2])

Science.spherical=spherical
Science.cartesian=cartesian


def rotation_matrix(ez,alpha):
    """
    Set a rotation matrix from the direction of the ez vector and a rotation angle alpha
    
    Parameter:
        ez: array (3)
            vector in the direction of the z-axis. 
            
        alpha: float (3) [rad]
            Rotation angle of the x-axis around z-axis (clockwise)
            
    Return:
        Msys2uni: float (3x3)
            Rotation matrix from the system defined by ez and the universal system.
            
        Muni2sys: float (3x3)
            Rotation matrix from the universal system to the system defined by ez
    """
    ez,one=spy.unorm(ez)
    ex=spy.ucrss([0,0,1],ez) #Spice is 5 faster for vcrss
    if spy.vnorm(ex)==0:
        ex=np.array([1,0,0]) if np.sum(ez)>0 else np.array([-1,0,0])
    ey=spy.ucrss(ez,ex)
    Msys2uni=np.array(list(np.vstack((ex,ey,ez)).transpose())).reshape((3,3))
    Muni2sys=spy.invert(Msys2uni)
    verbose("Rotation axis:",ex,ey,ez)
    return Msys2uni,Muni2sys

Science.rotation_matrix=rotation_matrix


