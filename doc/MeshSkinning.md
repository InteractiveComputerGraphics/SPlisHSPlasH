# MeshSkinning

MeshSkinning is a command line tool to generate a sequence of deformed meshes from a sequence of particle files of an elastic model. When simulating an elastic solid using SPlisHSPlasH, we only get particle data. If this data is exported using the PartioExported, the MeshSkinning tool is able to generate deformed triangle meshes in a post-processing step. The tool requires a triangle mesh of the reference configuration of the deformable solid. This mesh is then deformed according to the particle data. 

Below are two examples which were generated using the mesh skinning tool:

|![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/walrus_press.jpg)|![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/beam_meshes.jpg)|
| ---- | ---- |

The tool implements the methods of:

* Tassilo Kugelstadt, Jan Bender, José Antonio Fernández-Fernández, Stefan Rhys  Jeske, Fabian Löschner, and Andreas Longva. Fast Corotated Elastic SPH  Solids with Implicit Zero-Energy Mode Control. Proceedings of the ACM  on Computer Graphics and Interactive Techniques, 2021

## Parameters 

The mesh skinning tool can either use the scene file as input. Then you have to define the required triangle mesh in the scene file for each elastic model: 
```
"visMesh": "../models/beam.obj"
```

All transformations and required info is then automatically extracted from the scene file (assuming that the partio files are located in the output folder) which is the simplest way to use this tool.

Alternatively you can also define everything manually. Then you have to provide the partio path, mesh file as well as translation, rotation and scaling of the mesh file to fit the reference configuration of the particles. 

#### Command line options:

* -i, --input arg: Input file
* -o, --output arg: Output file
* -m, --mesh arg: Mesh file
* --scene arg: Scene file (all settings are imported from the scene file)
* --partioPath arg: Path of the partio files (when using a scene file). If not set, it is assumed that the files are in the standard output path.
* --scale arg: Scaling of input geometry (e.g. --scale "2 1 2") (default: 1 1 1)
* -t, --translation arg: Translation of input geometry (e.g. --translation "2 1 2") (default: 1 1 1)
* --axis arg: Rotation axis of input geometry (e.g. --axis "1 0 0") (default: 1 0 0)
* --angle arg: Angle of input geometry (e.g. --angle 1) (default: 0.0)
* -s, --startframe arg: Start frame (default: 1)
* -e, --endframe arg: End frame
* -r, --radius arg: Particle radius (default: 0.025)
* --supportRadiusFactor arg: The support radius is defined as factor*particleRadius (default: 6.0)
* --maxNeighbors arg: The maximum number of neighbors that are used for the interpolation. (default: 60)
* --splitting: Read a scene which used the object splitting export option.
* --overwrite: Overwrite existing files.
* -h, --help: Print help

#### Example: 
```
MeshSkinning --splitting --overwrite --scene ..\data\Scenes\Beam.json
```