# SPlisHSPlasH Scene Files

A SPlisHSPlasH scene file is a json file which can contain the following blocks:

* Configuration
* FluidBlocks
* FluidModels
* Emitters
* RigidBodies
* Fluid parameter block
* Animation fields

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
* cameraPosition (vec3): Initial position of the camera.
* cameraLookat (vec3): Lookat point of the camera.

##### Visualization:

* numberOfStepsPerRenderUpdate (int): Number of simulation steps per rendered frame
* renderWalls (int): 
  - 0: None
  - 1: Particles (all)
  - 2: Particles (no walls)
  - 3: Geometry (all)
  - 4: Geometry (no walls)

##### Export

* enablePartioExport (bool): Enable/disable partio export (default: false).
* enableVTKExport (bool): Enable/disable VTK export (default: false).
* enableRigidBodyExport (bool): Enable/disable rigid body export (default: false).
* dataExportFPS (float): Frame rate of particle and rigid body export (default: 25).
* particleAttributes (string): A list of attribute names separated by ";" that should be exported in the particle files (e.g. "velocity;density") (default: "velocity").
* enableStateExport (bool): Enable/disable export of complete simulation state (default: false).
* stateExportFPS (float): Frame rate of simulation state export (default: 1).

##### Simulation:

* timeStepSize (float): The initial time step size used for the time integration. If you use an adaptive time stepping, this size will change during the simulation (default: 0.001).
* particleRadius (float): The radius of the particls in the simulation (all have the same radius) (default: 0.025).
* sim2D (bool): If this parameter is set to true, a 2D simulation is performend instead of a 3D simulation (default: false).
* enableZSort (bool): Enable z-sort to improve cache hits and therefore to improve the performance (default: true).
* gravitation (vec3): Vector to define the gravitational acceleration (default: [0,-9.81,0]).
* maxIterations (int): Maximal number of iterations of the pressure solver (default: 100).
* maxError (float): Maximal density error in percent which the pressure solver tolerates (default: 0.01).
* boundaryHandlingMethod (int): The boundary handling method that is used in the simulation (default: 2, Volume Maps):
    - 0: particle-based boundaries (Akinci et al. 2012)
    - 1: density maps (Koschier et al. 2017)
    - 2: volume maps (Bender et al. 2019)
* simulationMethod (int): The pressure solver method used in the simulation (default: 4, DFSPH):
    - 0: Weakly compressible SPH for free surface flows (WCSPH)
    - 1: Predictive-corrective incompressible SPH (PCISPH)
    - 2: Position based fluids (PBF)
    - 3: Implicit incompressible SPH (IISPH)
    - 4: Divergence-free smoothed particle hydrodynamics (DFSPH)
    - 5: Projective Fluids (dynamic boundaries not supported yet)


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
        "rotationAxis": [0, 1, 0],
		"rotationAngle": 3.1415926535897932384626433832795,
        "velocity": 2,	
        "emitStartTime": 2,
		"emitEndTime": 6,
        "type": 0
    }
]
```

* type (int): Defines the shape of the emitter (default: 0).
  - 0: box
  - 1: circle
* width (int): Width of the box or radius of the circle emitter (default: 5).
* height (int): Height of the box (is only used for type 0) (default: 5).
* translation (vec3): Translation vector of the emitter (default: [0,0,0]).
* rotationAxis (vec3): Axis used to rotate the emitter. Note that in 2D simulations the axis is always set to [0,0,1] (default: [0,0,1]).
* rotationAngle (float): Rotation angle for the initial rotation of the emitter (default: 0). 
* velocity (float): Initial velocity of the emitted particles in direction of the emitter (default: 1).
* id: This id is used in the "Fluid parameter block" (see below) to define the properties of the fluid block. If no id is defined, then the standard id "Fluid" is used (default: "Fluid").	
* emitStartTime (float): Start time of the emitter (default: 0).
* emitEndTime (float): End time of the emitter (default: REAL_MAX).

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
        "isWall": true,
        "mapInvert": true, 
		"mapThickness": 0.0,
		"mapResolution": [20,20,20],
        "samplingMode": 1
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
* mapInvert (bool): Invert the map when using density or volume maps, flips inside/outside (default: false) 
* mapThickness (float): Additional thickness of a volume or density map (default: 0.0)
* mapResolution (vec3): Resolution of a volume or density map (defaut: [20,20,20])
* samplingMode (int): Surface sampling mode. 0 Poisson disk sampling, 1 Regular triangle sampling (default: 1).


## Fluid parameter block

```json
"Fluid":
{
    "density0": 1000, 
    "colorField": "velocity",
	"colorMapType": 1,
	"renderMinValue": 0.0,
	"renderMaxValue": 5.0,
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
    "emitterBoxMax": [0.0,4,4.0],
}
```

* density0 (float): Rest density of the corresponding fluid.

##### Particle Coloring 

* colorField (string): Choose vector or scalar field for particle coloring.
* colorMapType (int): Selection of a color map for coloring the scalar/vector field.
  - 0: None
  - 1: Jet
  - 2: Plasma
* renderMinValue (float): Minimal value used for color-coding the color field in the rendering process.
* renderMaxValue (float): Maximal value used for color-coding the color field in the rendering process.

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

##### Elasticity

* elasticityMethod (int): Elasticity method
  - 0: None
  - 1: Becker et al. 2009
  - 2: Peer et al. 2018
* youngsModulus (float): Young's modulus - coefficient for the stiffness of the material (default: 100000.0)
* poissonsRatio (float): Poisson's ratio - measure of the Poisson effect (default: 0.3)
* alpha (float): Coefficent for zero-energy modes suppression method (default: 0.0)
* elasticityMaxIter (float): (Peer et al. 2018) Maximum solver iterations (default: 100)
* elasticityMaxError (float): (Peer et al. 2019) Maximum elasticity error allowed by the solver (default: 1.0e-4)

##### Emitters

* maxEmitterParticles (int): Maximum number of particles the emitter generates. Note that reused particles (see below) are not counted here.
* emitterReuseParticles (bool):  Reuse particles if they are outside of the bounding box defined by emitterBoxMin, emitterBoxMax
* emitterBoxMin (vec3): Minimum coordinates of an axis-aligned box (used in combination with emitterReuseParticles)
* emitterBoxMax (vec3): Maximum coordinates of an axis-aligned box (used in combination with emitterReuseParticles)

## Animation fields

In this part the user can define one or more animation fields which animate fluid particles. The user can define math expressions for the components of the field quantity. The typical math terms like cos,sin,... can be used. 

Available expression variables:
- t: Current time.
- dt: Current time step size.
- x, y, z: Position of the particle which is in the animation field.
- vx, vy, vz: Velocity of the particle which is in the animation field.
- valuex, valuey, valuez: Value of the field quantity of the particle which is in the animation field.

Example:
```json
"particleField": "angular velocity",
"expression_x": "valuex + cos(2*t)"
```

This means that in each step we add cos(2*t) to the x-component of the angular velocity.

Example code:
```json
"AnimationFields": [
    {
        "particleField": "velocity",
        "translation": [-0.5, -0.5, 0],
        "rotationAxis": [0, 0, 1],
        "rotationAngle": 0.0,
        "scale": [0.5, 0.25, 0.8],
        "shapeType": 0,
        "expression_x": "cos(2*t)*0.1",
        "expression_y": "",
        "expression_z": ""
    }
]
```

* shapeType (int): Defines the shape of the animation field (default: 0).
  - 0: box
  - 1: sphere
  - 2: cylinder
* particleField (string): Defines the field quantity that should be modified by the field (e.g. velocity, angular velocity, position) (default: velocity)
* translation (vec3): Translation vector of the animation field (default: [0,0,0]).
* rotationAxis (vec3): Axis used to rotate the animation field (default: [0,0,1]).
* rotationAngle (float): Rotation angle for the initial rotation of the animation field (default: 0).
* scale (vec3): Scaling vector of the animation field.  
  - shapeType=0 (box): This vector defines the width, height, depth of the box.
  - shapeType=1 (sphere): The x-component of the vector defines the radius of the sphere. The other components are ignored.
  - shapeType=2 (cylinder): The x- and y-component of the vector defines the height and radius of the cylinder, repectively. The z-component is ignored.
* expression_x (string): Math expression for the x-component of the field quantity (default="").
* expression_y (string): Math expression for the y-component of the field quantity (default="").
* expression_z (string): Math expression for the z-component of the field quantity (default="").
