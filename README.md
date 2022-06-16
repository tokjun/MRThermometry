# MR Thermometry Module for 3D Slicer
Collection of 3D Slicer utility scripts for MR Thermometry


## Installation

The PRFThermometry module uses [scikit-image](https://scikit-image.org) for phase unwrapping. scikit-image
can be installed using pip from the Python interactor:

~~~~
import pip
pip.main(['install', 'scikit-image'])
~~~~


## Known issues
The color bar does not show up in the recent version of 3D Slicer due to the change in color bar management.

https://github.com/Slicer/Slicer/issues/4891
