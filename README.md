# Pryngles

## PlanetaRY spaNGLES

<!--[![PyPi version](https://pypip.in/v/pryngles/badge.png)](https://crate.io/packages/pryngles/)-->
<!--[![PyPi downloads](https://pypip.in/d/pryngles/badge.png)](https://crate.io/packages/pryngles/)-->
<!--Other badges: https://shields.io/category/activity -->

[![version](https://img.shields.io/pypi/v/pryngles?color=blue)](https://pypi.org/project/pryngles/)
[![downloads](https://img.shields.io/pypi/dw/pryngles)](https://pypi.org/project/pryngles/)
[![license](https://img.shields.io/pypi/l/pryngles)](https://pypi.org/project/pryngles/)
[![implementation](https://img.shields.io/pypi/implementation/pryngles)](https://pypi.org/project/pryngles/)
[![pythonver](https://img.shields.io/pypi/pyversions/pryngles)](https://pypi.org/project/pryngles/)
<!--[![codesize](https://img.shields.io/github/languages/repo-size/seap-udea/pryngles-public)](https://pypi.org/project/pryngles/)-->
[![arXiv](http://img.shields.io/badge/arXiv-2004.14121-orange.svg?style=flat)](http://arxiv.org/abs/2004.14121)

<!--
<p align="left">
<img src="https://github.com/seap-udea/pryngles-public/blob/master/gallery/pryngles-logo-wb.png?raw=true" alt="Logo" width="200"/>
</p>
-->

`Pryngles` is a `Python` package intended to produce useful
visualizations of the geometric configuration of a ringed exoplanet
(an exoplanet with a ring or exoring for short) and more importantly
to calculate the light-curve signatures produced by this kind of
planets.  The model behind the package has been developed in an effort
to predict the signatures that exorings may produce not only in the
light-curve of transiting exoplanets (a problem that has been
extensively studied) but also in the light of stars having
non-transiting exoplanets.

This is an example of what can be done with `Pryngles`:

<p align="center">
<img src="https://github.com/seap-udea/pryngles-public/blob/master/gallery/ecliptic-i_3.0e+01-lambobs_9.0e+01-betaobs_9.0e+01.gif?raw=true" alt="Animation" width="400"/>
</p>

For the science behind the model please refer to the following papers:

> Zuluaga, J.I., Sucerquia, M. & Alvarado-Montes, J.A. (2022), **The
  bright side of the light curve: a general photometric model for
  non-transiting exorings**, submitted to MNRAS (2022).

> Sucerquia, M., Alvarado-Montes, J. A., Zuluaga, J. I., Montesinos,
  M., & Bayo, A. (2020), **Scattered light may reveal the existence of
  ringed exoplanets**. Monthly Notices of the Royal Astronomical
  Society: Letters, 496(1), L85-L90.

<p align="center">
<img src="https://github.com/seap-udea/pryngles-public/blob/master/gallery/light-curve.png?raw=true" alt="Logo""/>
</p>

## Download and install

`pryngles` is available in `PyPI`, https://pypi.org/project/pryngles/.
To install it, just execute:

```
   pip install pryngles
```

If you prefer, you may download and install from the
[sources](https://pypi.org/project/pryngles/#files).

## Quick start

Any calculation in `Pryngles` start by creating a planetary system:

```python
import pryngles as pr
sys=pr.System()
```

Then we add objects to the planetary system using:

```python
S=sys.addStar()
P=sys.addPlanet(center=S,a=0.5,e=0.2,N=1000)
R=sys.addRing(center=P,fi=1.4,fe=2.4,i=30*pr.deg,N=1500)
O=sys.addObserver()
```

By default the planet has the radius of Saturn and it orbits a
solar-mass star.  In the example before the planet has a ring
extending from 1.4 to 2.4 planetary radius which is inclined 30
degrees with respect to the orbital plane. It has an orbit with
semimajor axis of 0.5 and eccentricity 0.2.  To discretize the surface
of the planet the program will use 1000 spangles or sequins (surface
elements with a circular shape) and the ring using 1500
spangles.

Once the system is set we can *ensamble* a simulation, ie. creating an
object able to produce a light-curve.

```python
RP=sys.ensambleSystem()
```

To see how the surface of the planet and the rings looks like run:

```python
RP.plotRingedPlanet()
```

You may change the position of the star in the orbit and see how the
appearance of the planet changes:

```python
RP.changeStellarPosition(45*pr.deg)
RP.plotRingedPlanet()
```

Below is the sequence of commands to produce your first light curve:

```python
import numpy as np
RP.changeObserver([90*pr.deg,30*pr.deg])
lambs=np.linspace(+0.0*pr.deg,+360*pr.deg,100)
Rps=[]
Rrs=[]
ts=[]
for lamb in lambs:
    RP.changeStellarPosition(lamb)
    ts+=[RP.t*sys.ut/pr.Const.days]
    RP.updateOpticalFactors()
    RP.updateDiffuseReflection()
    Rps+=[RP.Rip.sum()]
    Rrs+=[RP.Rir.sum()]

ts=np.array(ts)
Rps=np.array(Rps)
Rrs=np.array(Rrs)

#Plot
import matplotlib.pyplot as plt
fig=plt.figure()
ax=fig.gca()
ax.plot(ts,1e6*Rps,label="Planet")
ax.plot(ts,1e6*Rrs,label="Ring")
ax.plot(ts,1e6*(Rps+Rrs),label="Planet+Ring")

ax.set_xlabel("Time [days]")
ax.set_ylabel("Flux anomaly [ppm]")
ax.legend();
```

And *voilà*! 

<p align="center">
<img src="https://github.com/seap-udea/pryngles-public/blob/master/gallery/example-light-curve.png?raw=true" alt="Light curve"/>
</p>

Let's have some `Pryngles`.

## Tutorial

A complete `Jupyter` tutorial is provided
[here](https://github.com/seap-udea/pryngles-public/blob/master/pryngles-tutorial-exploration.ipynb).
If you prefer Google Colab check
[this](https://bit.ly/pryngles-tutorial-exploration).

In this tutorial you will learn how to configure, install, use and even modify `pryngles`.

## Disclaimer

<p align="center">
<img src="https://github.com/seap-udea/pryngles-public/blob/master/gallery/disco-planet.jpeg?raw=true" alt="Logo" width="150"/>
</p>

This is the *disco* version of Pryngles.  We are improving the
resolution and performance of the software for future releases.

## What's new

- **0.4.x versions**:

  - A new model to create and populate planetary system has been
    implemented.

- **0.3.x versions**:

  - A water mark with version number included.
  - Version is now available in the __version__ variable.
  - Scattering formulae tested and verified.
  - Package has been compared against similar packages (good
    agreement) but disclaimer has been done.
  - New version number scheme: 0.x.y (x-major, y-minor release),
    0.x.y.z (z test version).
  - Major corrections in diffuse formulae.

- **0.2.1.x versions**:

  - Tutorial is now working in Google Colab.
  - References were corrected.
  - The home url was set as the PyPI web page.
  - Non-linear (4th order) limb darkening included.
  - Added the class `Extra`.
  - Function to draw logo: `drawPryngles`.
  - Added function `prynglesMark`.
  - Now `__version__` variable is available.

- **0.2.0.x versions**:

  - First official version of the package.

------------

This package has been designed and written originally by Jorge
I. Zuluaga, Mario Sucerquia & Jaime A. Alvarado-Montes (C) 2022
