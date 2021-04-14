# Paraview Partio Plugin

The plugin can be used to import [Partio](https://www.disneyanimation.com/technology/partio.html)  particle data in Paraview. This plugin can be used to import and render the particle data generated with SPlisHSPlasH.

## Installation

1. This add-on requires the partio python module. To build this module first install pybind11 by 
	
   	pip install pybind11

   Then the module can be built by the calling

   	python setup.py build_ext

   in the directory partio_extension. Note that you have to use the same Python version as your Paraview version uses (for Paraview 5.9.0 this is Python 3.8).

3. Copy the generated file partio.* (name depends on system) to the Paraview site-packages folder (in Windows: bin\Lib\site-packages).

4. Start Paraview and load plugin (Tools/Manage Plugins.../Load New...)

4. Optinally activate Auto Load.

## Usage

After loading the plugin you can directly open partio files.
