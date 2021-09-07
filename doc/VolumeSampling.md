# VolumeSampling

The simulator can load particle data from partio files. This particle data then defines the initial configuration of the particles in the simulation. The VolumeSampling tool allows you to sample a volumetric object with particle data. This means you can load an OBJ file with a closed surface geometry and sample the interior with particles using different methods. Especially when simulating elastic solids a good sampling is beneficial as shown by Kugelstadt et al. [2021].

Below are two examples which were generated using the volume sampling tool:

| ![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/volume_sampling1.jpg) | ![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/volume_sampling2.jpg) |
| ------------------------------------------------------------ | ------------------------------------------------------------ |

The tool implements the methods of:

* M. Jiang, Y. Zhou, R. Wang, R. Southern, J. J. Zhang. Blue noise sampling using an SPH-based method. ACM Transactions on Graphics, 2015
* Tassilo Kugelstadt, Jan Bender, José Antonio Fernández-Fernández, Stefan Rhys  Jeske, Fabian Löschner, and Andreas Longva. Fast Corotated Elastic SPH  Solids with Implicit Zero-Energy Mode Control. Proceedings of the ACM  on Computer Graphics and Interactive Techniques, 2021

#### Command line options:

* -h, --help: Print help
* -i, --input arg: Input file (obj)
* -o, --output arg: Output file (bgeo or vtk)
* -r, --radius arg: Particle radius (default: 0.025)
* -s, --scale arg: Scaling of input geometry (e.g. --scale "2 1 2") (default: 1 1 1)
* -m, --mode arg: Mode (regular=0, almost dense=1, dense=2, Jiang et al. 2015=3, Kugelstadt et al. 2021=4) (default: 4)
* --region arg: Region to fill with particles (e.g. --region "0 0 0 1 1 1")
* --steps arg: SPH time steps (default: 100)
* --cflFactor arg: CFL factor (default: 0.25)
* --viscosity arg: Viscosity coefficient (XSPH) (default: 0.25)
* --cohesion arg: Cohesion coefficient
* --adhesion arg: Adhesion coefficient
* --stiffness arg: Stiffness coefficient (only mode 3) (default: 10000.0)
* --dt arg: Time step size (only mode 3) (default: 0.0005)
* --res arg: Resolution of the Signed Distance Field (e.g. --res "30 30 30")
* --invert: Invert the SDF to sample the outside of the object in the bounding box/region
* --no-cache: Disable caching of SDF.

#### Example: 
```
VolumeSampling.exe --mode 4 -i ..\data\models\bunny.obj -o bunny.vtk
```

