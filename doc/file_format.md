# SPlisHSPlasH Scene Files

A SPlisHSPlasH scene file is a json file which can contain the following blocks:

* Configuration
* FluidBlocks
* FluidModels
* Emitters
* RigidBodies
* Fluid parameter block

## Configuration

This part contains the general settings of the simulation and the pressure solver. 

Example code:
```json
"Configuration": 
{
    "pause": true,
	"sim2D": false, 
    "timeStepSize": 0.001,
    "numberOfStepsPerRenderUpdate": 2,
    "colorMapType": 1,
    "particleRadius": 0.025, 
    "simulationMethod": 4,
    "gravitation": [0.0,-9.81,0], 
    "cflMethod": 1, 
    "cflFactor": 1,
    "cflMaxTimeStepSize": 0.005,
    "maxIterations": 100,
    "maxError": 0.01,
    "maxIterationsV": 100,
    "maxErrorV": 0.1,		
    "stiffness": 50000,
    "exponent": 7,
    "velocityUpdateMethod": 0,
    "enableDivergenceSolver": true
}
```
##### General:

* pause (bool): Pause simulation at beginning.
* pauseAt (float): Pause simulation at the given time. When the value is negative, the simulation is not paused.
* stopAt (float): Stop simulation at the given time and exit. When the value is negative, the simulation is not stopped.

##### Visualization:

* numberOfStepsPerRenderUpdate (int): Number of simulation steps per rendered frame
* colorField (int): Choose vector or scalar field for particle coloring.
* colorMapType (int): Selection of a color map for coloring the scalar/vector field.
  - 0: None
  - 1: Jet
  - 2: Plasma
* renderMinValue (float): Minimal value used for color-coding the color field in the rendering process.
* renderMaxValue (float): Maximal value used for color-coding the color field in the rendering process.
* renderWalls (int): 
  - 0: None
  - 1: Particles (all)
  - 2: Particles (no walls)
  - 3: Geometry (all)
  - 4: Geometry (no walls)

##### Export

* enablePartioExport (bool): Enable/disable partio export.
* partioFPS (int): Frame rate of partio export.

##### Simulation:

* timeStepSize (float): The initial time step size used for the time integration. If you use an adaptive time stepping. This size will change during the simulation.
* particleRadius (float): The radius of the particls in the simulation (all have the same radius).
* sim2D (bool): If this parameter is set to true, a 2D simulation is performend instead of a 3D simulation.
* gravitation (vec3): Vector to define the gravitational acceleration.
* maxIterations (int): Maximal number of iterations of the pressure solver.
* maxError (float): Maximal density error in percent which the pressure solver tolerates.
* simulationMethod (int): The pressure solver method used in the simulation:
    - 0: Weakly compressible SPH for free surface flows (WCSPH)
    - 1: Predictive-corrective incompressible SPH (PCISPH)
    - 2: Position based fluids (PBF)
    - 3: Implicit incompressible SPH (IISPH)
    - 4: Divergence-free smoothed particle hydrodynamics (DFSPH)
    - 5: Projective Fluids


##### WCSPH parameters:

* stiffness (float): Stiffness coefficient of the equation of state.
* exponent (float): Exponent in the equation of state.

##### PBF parameters:

* velocityUpdateMethod (int): 
  - 0: First Order Update
  - 1: Second Order Update

##### DFSPH parameters:

* enableDivergenceSolver (bool): Turn divergence solver on/off.
* maxIterationsV (int): Maximal number of iterations of the divergence solver.
* maxErrorV (float): Maximal divergence error in percent which the pressure solver tolerates.

##### Projective Fluids parameters:

* stiffness (float): Stiffness coefficient used by the pressure solver.

##### Kernel: 

* kernel (int): Kernel function used in the SPH model. 
  - For a 3D simulation:
    - 0: Cubic spline
    - 1: Wendland quintic C2
    - 2: Poly6
    - 3: Spiky
    - 4: Precomputed cubic spline (faster than cubic spline)
  - For a 2D simulation:
    - 0: Cubic spline
    - 1: Wendland quintic C2
* gradKernel (int): Gradient of the kernel function used in the SPH model.
  - For a 3D simulation:
    - 0: Cubic spline
    - 1: Wendland quintic C2
    - 2: Poly6
    - 3: Spiky
    - 4: Precomputed cubic spline (faster than cubic spline)
  - For a 2D simulation:
    - 0: Cubic spline
    - 1: Wendland quintic C2


##### CFL:

* cflMethod (int): CFL method used for adaptive time stepping.  
  - 0: No adaptive time stepping 
  - 1: Use CFL condition
  - 2: Use CFL condition and consider number of pressure solver iterations
* cflFactor (float): Factor to scale the CFL time step size.
* cflMaxTimeStepSize (float): Max. allowed time step size.

## FluidBlocks

In this part the user can define multiple axis-aligned blocks of fluid particles.

Example code:
```json
"FluidBlocks": [
    {
        "denseMode": 0,
        "start": [-2.0, 0.0, -1],
        "end": [-0.5, 1.5, 1],
        "translation": [1.0, 0.0, 0.0],
        "scale": [1, 1, 1]
    }
]	
```

* start (vec3): Minimum corrdinate of the box which defines the fluid block.
* end (vec3): Maximum corrdinate of the box which defines the fluid block.
* translation (vec3): Translation vector of the block.
* scale (vec3): Scaling vector of the block. 
* denseMode (int): 
  - 0: regular sampling
  - 1: more dense sampling
  - 2: dense sampling
* initialVelocity (vec3): The initial velocity is set for all particles in the block.
* id: This id is used in the "Fluid parameter block" (see below) to define the properties of the fluid block. If no id is defined, then the standard id "Fluid" is used.

## FluidModels

This part can be used to import one or more partio particle files in the scene.

Example code:
```json
"FluidModels": [
    {
        "particleFile": "../models/bunny.bgeo",
        "translation": [-2.0, 0.1, 0.0],
        "rotationAxis": [0, 1, 0],
        "rotationAngle": 3.14159265359,
        "scale": 1
    }
]
```

* particleFile (string): Path of the partio file which contains the particle data.
* translation (vec3): Translation vector of the fluid model.
* scale (vec3): Scaling vector of the fluid model. 
* rotationAxis (vec3): Axis used to rotate the particle data after loading.
* rotationAngle (float): Rotation angle for the initial rotation of the particle data. 
* id: This id is used in the "Fluid parameter block" (see below) to define the properties of the fluid block. If no id is defined, then the standard id "Fluid" is used.

## Emitters

In this part the user can define one or more emitters which generate fluid particles.

Example code:
```json
"Emitters": [
    {
        "width": 5, 
        "height": 5, 
        "translation": [-1,0.75,0.0],
        "direction": [1, 0, 0],
        "velocity": [2, 0, 0.0],	
        "comment" : "The initial velocity should be at least 2*particleRadius*emitsPerSecond.",
        "emitsPerSecond": 40,
        "type": 0
    }
]
```

* type (int): Defines the shape of the emitter.
  - 0: box
  - 1: circle
* width (int): Width of the box or radius of the circle emitter.
* height (int): Height of the box (is only used for type 0).
* translation (vec3): Translation vector of the emitter.
* direction (vec3): Direction of the emitter. This defines the normal vector of the box/circle.
* velocity (vec3): Initial velocity of the emitted particles.
* emitsPerSecond (float): Defines the number of emits per second.
* id: This id is used in the "Fluid parameter block" (see below) to define the properties of the fluid block. If no id is defined, then the standard id "Fluid" is used.	

## RigidBodies

Here, the static and dynamic rigid bodies are defined which define the boundary in the scene. 
In case of dynamic rigid bodies, the PositionBasedDynamics library is used for their simulation. Note that in this case the PositionBasedDynamics library also reads this json scene files and picks out the relevant parts. That means if you want to define for example a hinge joint or a motor, then just use the json format of PositionBasedDynamics in this scene file.

Example code:
```json
"RigidBodies": [
    {
        "geometryFile": "../models/UnitBox.obj",
        "translation": [0,2,0],
        "rotationAxis": [1, 0, 0],
        "rotationAngle": 0,
        "scale": [2.5, 4, 1.0],
        "color": [0.1, 0.4, 0.6, 1.0], 
        "isDynamic": false,
        "isWall": true
    }
]
```

* geometryFile (string): Path to a OBJ file which contains the geometry of the body.
* particleFile (string): Path to a partio file which contains a surface sampling of the body. Note that the surface sampling is done automatically if this parameter is missing.
* translation (vec3): Translation vector of the rigid body.
* scale (vec3): Scaling vector of the rigid body.
* rotationAxis (vec3): Axis used to rotate the rigid body after loading.
* rotationAngle (float): Rotation angle for the initial rotation of the rigid body. 
* isDynamic (bool): Defines if the body is static or dynamic.
* isWall (bool): Defines if this is a wall. Walls are typically not rendered. This is the only difference.
* color (vec4): RGBA color of the body.



## Fluid parameter block

```json
"Fluid":
{
    "density0": 1000, 
    "surfaceTension": 0.2,
    "surfaceTensionMethod": 0,		
    "viscosity": 0.01,
    "viscosityMethod": 1, 
    "vorticityMethod": 1, 
    "vorticity": 0.15, 
    "viscosityOmega": 0.05,
    "inertiaInverse": 0.5,
    "maxEmitterParticles": 1000,
	"emitterReuseParticles": false,
	"emitterBoxMin": [-4.0,-1.0,-4.0],
	"emitterBoxMax": [0.0,4,4.0]
}
```

* density0 (float): Rest density of the corresponding fluid.

##### Viscosity

* viscosityMethod (int): Viscosity method
  - 0: None
  - 1: Standard
  - 2: XSPH
  - 3: Bender and Koschier 2017
  - 4: Peer et al. 2015
  - 5: Peer et al. 2016
  - 6: Takahashi et al. 2015 (improved)
  - 7: Weiler et al. 2018
* viscosity (float): Coefficient for the viscosity force computation
* viscoMaxIter (int): (Implicit solvers) Max. iterations of the viscosity solver.
* viscoMaxError (float): (Implicit solvers) Max. error of the viscosity solver.
* viscoMaxIterOmega (int): (Peer et al. 2016) Max. iterations of the vorticity diffusion solver.
* viscoMaxErrorOmega (float): (Peer et al. 2016) Max. error of the vorticity diffusion solver.
* viscosityBoundary (float): (Weiler et al. 2018) Coefficient for the viscosity force computation at the boundary.


##### Vorticity

* vorticityMethod (int): Vorticity method
  - 0: None
  - 1: Micropolar model
  - 2: Vorticity confinement
* vorticity (float): Coefficient for the vorticity force computation
* viscosityOmega (float): (Micropolar model) Viscosity coefficient for the angular velocity field.
* inertiaInverse (float): (Micropolar model) Inverse microinertia used in the micropolar model.

##### Drag force

* dragMethod (int): Drag force method
  - 0: None
  - 1: Macklin et al. 2014
  - 2: Gissler et al. 2017
* drag (float): Coefficient for the drag force computation

##### Surface tension

* surfaceTensionMethod (int): Surface tension method
  - 0: None
  - 1: Becker & Teschner 2007
  - 2: Akinci et al. 2013
  - 3: He et al. 2014
* surfaceTension (float): Coefficient for the surface tension computation