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
<!--[![arXiv](http://img.shields.io/badge/arXiv-2004.14121-orange.svg?style=flat)](http://arxiv.org/abs/2004.14121)-->
[![arXiv](http://img.shields.io/badge/arXiv-2207.08636-orange.svg?style=flat)](http://arxiv.org/abs/2207.08636)
[![ascl](https://img.shields.io/badge/ascl-2205.016-blue.svg?colorB=262255)](https://ascl.net/2205.016)

<!--
<p align="left">
<img src="https://github.com/seap-udea/pryngles-public/blob/master/gallery/pryngles-logo-wb.png?raw=true" alt="Logo" width="200"/>
</p>
-->

`Pryngles` is a `Python` package intended to produce useful
visualizations of the geometric configuration of a ringed exoplanet
(an exoplanet with a ring or exoring for short) and more importantly
to calculate the light curve produced by this kind of planets.  The
model behind the package has been developed in an effort to predict
the signatures that exorings may produce not only in the light curve
of transiting exoplanets (a problem that has been extensively studied)
but also in the light of stars having non-transiting exoplanets (the
bright side of the light curve).

*If `PyPI` does not render properly the images if this README please
check it in our [public github
repo](https://github.com/seap-udea/pryngles-public).*

This is an example of what can be done with `Pryngles`:

<p align="center">
<img src="https://raw.githubusercontent.com/seap-udea/pryngles-public/master/gallery/ecliptic-i_3.0e%2B01-lambobs_9.0e%2B01-betaobs_9.0e%2B01.gif" alt="Animation" width="400"/>
</p>

For the science behind the model please refer to the following papers:

> Zuluaga, J.I., Sucerquia, M. & Alvarado-Montes, J.A. (2022), **The
  bright side of the light curve: a general photometric model for
  non-transiting exorings**, accepted for publication in Astronomy and
  Computing (2022), [arXiv:2207.08636](https://arxiv.org/abs/2207.08636).

> Sucerquia, M., Alvarado-Montes, J. A., Zuluaga, J. I., Montesinos,
  M., & Bayo, A. (2020), **Scattered light may reveal the existence of
  ringed exoplanets**. Monthly Notices of the Royal Astronomical
  Society: Letters, 496(1), L85-L90.

<p align="center">
<img src="https://raw.githubusercontent.com/seap-udea/pryngles-public/master/gallery/light-curve.png" alt="Logo""/>
</p>

## Download and install

`pryngles` is available in `PyPI`, https://pypi.org/project/pryngles/.
To install it, just execute:

```
   pip install -U pryngles
```

If you prefer, you may download and install from the
[sources](https://pypi.org/project/pryngles/#files).

## Quick start

Import the package and some useful utilities:

```python
import pryngles as pr
from pryngles import Consts
```

> **NOTE**: If you are working in `Google Colab` before producing any plot please load the 
  matplotlib backend:

  ```python
  %matplotlib inline
  ```

Any calculation in `Pryngles` starts by creating a planetary system:

```python
sys=pr.System()
```

Then we add objects to the planetary system using:

```python
S=sys.add(kind="Star",
          physics=dict(radius=Consts.rsun/sys.ul),
          optics=dict(limb_coeffs=[0.65])
         )
P=sys.add(kind="Planet",primary=S,
          orbit=dict(a=0.2,e=0.0),
          physics=dict(radius=Consts.rsaturn/sys.ul)
         )
R=sys.add(kind="Ring",primary=P,
          physics=dict(fi=1.5,fe=2.5,i=30*Consts.deg)
         )
O=sys.add(kind="Observer",
          optics=dict(lamb=90*Consts.deg,beta=90*Consts.deg)
         )
```

In the example before the planet has a ring extending from 1.5 to 2.5
planetary radius which is inclined 30 degrees with respect to the
orbital plane. It has an orbit with semimajor axis of 0.2 and
eccentricity 0.0.

Once the system is set we can *ensamble* a simulation, ie. creating an
object able to produce a light-curve.

```python
RP=sys.ensamble_system()
```

To see how the surface of the planet and the rings looks like run:

```python
RP.plotRingedPlanet()
```

You may change the position of the star in the orbit and see how the
appearance of the planet changes:

```python
RP.changeStellarPosition(45*Consts.deg)
RP.plotRingedPlanet()
```

Below is the sequence of commands to produce your first light curve:

```python
import numpy as np
RP.changeObserver([90*Consts.deg,30*Consts.deg])
lambs=np.linspace(+0.0*Consts.deg,+360*Consts.deg,100)
Rps=[]
Rrs=[]
ts=[]
for lamb in lambs:
    RP.changeStellarPosition(lamb)
    ts+=[RP.t*sys.ut/Consts.day]
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
ax.plot(ts,Consts.ppm*Rps,label="Planet")
ax.plot(ts,Consts.ppm*Rrs,label="Ring")
ax.plot(ts,Consts.ppm*(Rps+Rrs),label="Planet+Ring")

ax.set_xlabel("Time [days]")
ax.set_ylabel("Flux anomaly [ppm]")
ax.legend();
```

And *voilà*! 

<p align="center">
<img src="https://raw.githubusercontent.com/seap-udea/pryngles-public/master/gallery/example-light-curve.png" alt="Light curve"/>
</p>

Let's have some `Pryngles`.

## Tutorials

We have prepared several `Jupyter` tutorials to guide you in the usage
of the package. The tutorials evolve as the package is being optimized.

- **Quickstart**.  In this tutorial you will learn the basics about the package. 
  [Download](https://github.com/seap-udea/pryngles-public/blob/master/pryngles-tutorial-quickstart.ipynb), [Google Colab](https://drive.google.com/file/d/12ZzuhVlKNaxJcKhGAPKScHvkWQYa_fAJ/view?usp=sharinghttps://bit.ly/pryngles-tutorials-quickstart).

## Disclaimer

<p align="center">
<img src="https://raw.githubusercontent.com/seap-udea/pryngles-public/master/gallery/disco-planet.jpeg" alt="Logo" width="150"/>
</p>

This is the *disco* version of Pryngles.  We are improving resolution,
performance, modularity and programming standards for future releases.

## What's new
   
- **0.6.x versions**:

  - 0.6.0 is the official release version, after paper acceptance and
    arXiv submission.
  - Link to quickstart tutorial in Google Colab, updated.
  - Updated information about paper in the arXiv and ACL code.

- **0.5.x versions**:
  
  - Preview method plotRingedPlanet modified to work under Google Colab.
  - Physical and astronomical constants included.
  - A new tutorial was included.
  - A major update in the classes to create and populate planetary
    system.

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
