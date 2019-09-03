# SPlisHSPlasH - Maya Plugin

## Installation

To install the plugin, simply add the path of the plugin (SPlisHSPlasH.py) to the environment variables MAYA_PLUG_IN_PATH and MAYA_SCRIPT_PATH. 

The plugin can be loaded in the "Plug-In Managar" and then a new menu "SPlisHSPlasH" appears.

## Generating a scene

#### Initializing a scene

* Add at least one scene configuration.
* Add at least one fluid material.

#### Convert selection to fluid/rigid bodies

* Import or generate geometry that should be used as fluid boundary or rigid body.
* Select the corresponding geometry nodes.
* Click on "Convert selection to fluid" or "Convert selection to rigid bodies" in the menu.
* Select the transform of such a geometry to adapt the settings (in "Extra Attributes").

**Important**: When generating a fluid, the geometry should be a closed mesh since the simulator will sample the interior of the geometry.

#### Create emitter

* Click on "Create rectangular emitter" or "Create circular emitter".
* Adapt the transformation as needed. 
* The scaling and the velocity of the emitter defines how many particles are emitted per step.

**Important**: The emitter geometry should not be modified (except scaling). The geometry is just used to give a visual feedback and not to generate the scene.

#### Create animation field

* Create an animation field by clicking on the corresponding menu item. 
* Adapt the transformation as needed.
* Add an expression in the "Extra Attributes".

**Important**: The emitter geometry should not be modified (except scaling). The geometry is just used to give a visual feedback and not to generate the scene.

#### Save the scene

* Click on "Save scene".
* The scene is stored as json file and all used geometries are exported in the same directory as obj files.
* The scene can be directly used by the SPlisHSPlasH simulator.

## Importing rigid body data

#### Generating rigid body data in a simulation

When running a simulation in SPlisHSPlasH you can export the rigid body geometries and motions by setting the flag "enableRigidBodyExport" to true. This can also be done during runtime. The data is stored in the output directoy.

#### Importing rigid body data in Maya

First, load the SPlisHSPlasH plugin in the "Plug-In Managar". In the SPlisHSPlasH menu select "Load rigid body data". A file dialog appears. In this dialog you should select the first file of your exported rigid body data. The data is then imported completely.

**Important**: Always select the first file of a sequence since it contains additional data. 

