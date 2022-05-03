# Pryngles

## PlanetaRY spaNGLES

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
  non-transiting exorings**, in preparation (2022).

> Sucerquia, M., Alvarado-Montes, J. A., Zuluaga, J. I., Montesinos,
  M., & Bayo, A. (2020), **Scattered light may reveal the existence of
  ringed exoplanets**. Monthly Notices of the Royal Astronomical
  Society: Letters, 496(1), L85-L90.

## Download and install

`pryngles` is available in `PyPI`, https://pypi.org/project/pryngles/.
To install it, just execute:

```
   pip install pryngles
```

If you prefer, you may download and install from the
[sources](https://pypi.org/project/pryngles/#files).

## Tutorial

A complete `Jupyter` tutorial is provided
[here](https://github.com/seap-udea/pryngles-public/blob/master/pryngles-tutorial-exploration.ipynb).
If you prefer Google Colab check
[this](https://bit.ly/pryngles-tutorial-exploration).

In this tutorial you will learn how to configure, install, use and even modify `pryngles`.

## What's new

- **0.2.0.x versions**:

  - First official version of the package.

- **0.2.1.x versions**:

  - Tutorial is now working in Google Colab.
  - References were corrected.
  - The home url was set as the PyPI web page.
  - Non-linear (4th order) limb darkening included.
  - Added the class `Extra`.
  - Function to draw logo: `drawPryngles`.
  - Added function `prynglesMark`.
  - Now `__version__` variable is available.

------------

This package has been designed and written originally by Jorge
I. Zuluaga, Mario Sucerquia & Jaime A. Alvarado-Montes (C) 2022
