# Software Architecture

![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/SoftwareArchitechture.jpg)

SPlisHSPlasH follows a very intuitive and modular design approach. We want to illustrate part of the software architecture in conjunction with the simplified class diagram above. Note, that this documentation only covers the simulation part of SPlisHSPlasH. The whole software architecture follows a similar design pattern as the **Model View Controller**. 

## The Simulation class

The simulation class is the main part of the software. It contains the currently used simulation method (`TimeStep`), all fluids (`FluidModel`), all boundaries (`BoundaryModel`), and a `AnimationFieldSystem`. It is defined as a **singleton**, thus only one simulation instance exists during the runtime. The simulation instance contains:
- exactly one `TimeStep` instance, which defines the simulation loop and contains the pressure solver
- any number of `FluidModel` instances each defining a different fluid phase
- any number of `BoundaryModel` instances representing either dynamic rigid bodies or static boundaries 
- exactly one `AnimationFieldSystem` instance which allows to animate particles in a predefined area

The simulation class also implements the following:
- evaluation of the SPH kernel methods
- update of the time step size using a CFL condition
- uniform invocation of all `EmitterSystem` instances
- invocation of `AnimationFieldSystem` instance
- saving & loading the current simulation state

Lastly, the simulation class also contains a well defined interface for the neighborhood search functionalities defined in [CompactNSearch](https://github.com/InteractiveComputerGraphics/CompactNSearch) or [cuNSearch](https://github.com/InteractiveComputerGraphics/cuNSearch), which are further needed in the respective algorithm implementations in e.g. the TimeStep or NonPressureForces.

## The TimeStep class

The TimeStep class is a **abstract base class** for any subsequent derived simulation method one wants to implement. It implements the required interface for the simulation class, noteably the `step()` function containing the simulation algorithm called in the main loop. During execution there exists exactly one instance of a TimeStep class. By default SPlisHSPlasH currently implements the following pressure solvers and the corresponding simulation algorithms:

- WCSPH
- PCISPH
- PBF
- IISPH
- DFSPH
- Projective Fluids

## The FluidModel class

A FluidModel instance represents a fluid phase with its respective properties and applied effects to it. SPlisHSPlasH allows for arbitrary many FluidModels inside a simulation as long as there is at least one and they all have a different id (see scene file format). One FluidModel contains the following:

- Physical parameters like rest density, mass, position, velocity, acceleration and current density
- Simulation parameters like the number of particles, their state and ID
- References to the applied non-pressure effects, one for each: 
  - Drag
  - Elasticity
  - Surface tension
  - Viscosity
  - Vorticity
- Emitter systems 

Concerning the non-pressure effects, each FluidModel can *only* utilize up to one method per non-pressure effect, which will be directly included in the computation inside the `computeNonPressureForces()` method of the `Simulation` class. Thus having e.g. two different surface tension algorithms inside one FluidModel is **not** possible. However, it is possible to define e.g. two phases, which have a different viscosity model and only one regarding surface tension effects. 

The emitters are only stored inside the FluidModels since they are assigned to a fixed FluidModel. Their functionalities are uniformly executed by the Simulation class in the `emitParticles()` step usually invoked at the end of the simulation loop of the current TimeStep instance.

## The BoundaryModel class

The BoundaryModel class provides a useful base class for any boundary handling methods. It stores a `RigidBodyObject` reference representing the object of the boundary. This can be a stationary or dynamic rigid body, whose coupling effects are handled uniformly. Note that `RigidBodyObject` is an abstract class providing an interface for the two derived classes `StaticRigidBody` and `PBDRigidBody`. The first is handled internally and represent stationary objects. The latter describes a moving rigid body which is simulated externally by the [PositionBasedDynamics](https://github.com/InteractiveComputerGraphics/PositionBasedDynamics) library. SPlisHSPlasH implements three different boundary models:

- Particle-based rigid-fluid coupling [Akinci et al. 2012]
- Density maps [Koschier and Bender 2017]
- Volume maps [Bender et al. 2019]

Finally, SPlisHSPlasH defines a boundary as a list of rigid bodies in conjunction with a rigid-fluid coupling algorithm.

