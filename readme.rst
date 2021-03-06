GAMS - Geoscience Analyst Magnetics Suite
====

GAMS is a lightweight GUI to read, visualize, and operate on magnetic data imported from Geoscience Analyst. 
The GUI is designed to import magnetic data from a Geoscience Analyst (GA) workspace, perform operations and transformations on the magnetic grids (e.g., spatial derivates, tilt transform), and export the results back into GA for final visualization. Grid pre-processing options are given, and fast (and lazy) vectorized calculations allow visualization and selection of optimal settings in real time.

Installation
====

GAMS can be installed with pip + git using::

	pip install git+https://github.com/eroots/gams.git

Alternatively, the repository can be manually downloaded and installed using the install script, i.e., by navigating to the GAMS folder and running::

	python setup.py install

Geoscience Analyst (a free geoscientific data viewer from Mira Geoscience) can be installed from:

https://mirageoscience.com/mining-industry-software/geoscience-analyst/


Dependencies
====

* matplotlib
* numpy
* scipy
* pyqt
* geoh5py

Installing GAMS through pip + git will result in these libraries being installed in your current python environment.

Getting Started
====

Once installed, the GUI can be launched from the command line as::

	gams

A file dialog box will open where you can load the Geoscience Analyst workspace (.geoh5 file) you wish to work from. A second dialog box will open where you can select the specific grid you want to work on. Note that currently all grids from the selected workspace will be displayed and there is no check that the grid you select is magnetic data.

Test data is available in the gams/resources folder.

Documentation
====

Documentation, including explanations of the pre-processing options, is included in gams.pdf