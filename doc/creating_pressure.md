# Creating Pressure Solvers

SPlisHSPlasH organizes the pressure solvers in their respective folders inside the `/SPlisHSPlasH/` directory. For example DFSPH can be found inside `/SPlisHSPlasH/DFSPH/`. We highly suggest the user to follow our file organization scheme. The user can also add new pressure solvers by by creating new or copying and modifying existing classes and then adding them to the build system plus additionally registering in the source code. 

Note that we do not strictly distinguish the pressure solver from the simulation algorithm. Each `TimeStep` class implements a whole time step including the pressure solver. The non-pressure forces are decoupled in their respective classes and only implicitly called. Thus for implementing a new pressure solver, we suggest copying the files from for example WCSPH and replacing the pressure solver by your own one. Note further, that we usually decouple data from the algorithm with the `SimulationData` classes. We strongly recommend doing the same with your implementation.

## Creating a new class

Again, we want to stress that copying and modifying existent methods is easier than writing a new class from scratch. However, if you want to do so, be sure to implement every abstract method inherited from `TimeStep`. These include:
- `void step()`, the simulation step function
- `void resize()`,  a method to initialize and resize any used field 

Albeit being not necessary, the user may also want to override/redefine the following methods:
- `void init()`, the initialization method. It is **important to call** `TimeStep::init()` inside this method
- `void reset()`, the method invoked on every reset command
- `void computeDensities()`, if the user does not want to utilize the given density computation 

A minimal working example of a derive class is shown below:

**TimeStepMyPressureSolver.h**

```cpp
#ifndef __TimeStepMyPressureSolver_h__
#define __TimeStepMyPressureSolver_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class TimeStepMyPressureSolver : public TimeStep
	{
    public:
        TimeStepMyPressureSolver();
        virtual ~TimeStepMyPressureSolver();

        virtual void step();
        
        virtual void resize();
    };
}

#endif
```

**TimeStepMyPressureSolve.cpp**

```cpp
#include "TimeStepMyPressureSolve.h"

using namespace SPH;
using namespace GenParam;

TimeStepMyPressureSolve::TimeStepMyPressureSolve() :
    TimeStep()
{
    [...]
}

TimeStepMyPressureSolve::~TimeStepMyPressureSolve(void)
{
    [...]
}

void TimeStepMyPressureSolve::step()
{
    [...]
}

void TimeStepMyPressureSolve::resize()
{
    [...]
}
```

SPlisHSPlasH assumes your simulation method allows for operator splitting, thus usually dividing the simulation into non-pressure forces and the pressure solver plus advection. The latter is subject of the TimeStep class. It is still possible to implement these together inside your own TimeStep class, but it contradicts SPlisHSPlasH's design principles. Since the `step()` method is forwarded to the main loop by the simulation class, its purpose is to define the simulation algorithm. For guidance, we also provide a simple SPH simulation algorithm outline:

```cpp
void TimeStepWCSPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

    // 1. Perform a neighborhood search
	performNeighborhoodSearch();

    // 2. Compute non-pressure forces and SPH densities
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		clearAccelerations(fluidModelIndex);
		computeDensities(fluidModelIndex);
	}
	sim->computeNonPressureForces();

    // 3. Compute pressure forces
	computePressureForces();

    // 4. Update time step tize with CFL condition
	sim->updateTimeStepSize();

    // 5. Advect particles
	advectParticles();

    // 6. Emit and/or animate particles if necessary
	sim->emitParticles();
	sim->animateParticles();

    // 7. Advect time
	tm->setTime(tm->getTime() + h);
}
```
where `computeDensities(...)` and `clearAcceleration(...)` are already defined by the base class.

We recommend the user to split the simulation algorithm and its data into two separate classes as it is the case for our already implemented ones. 

## Registering the pressure solver

To add our new simulation method, we have to integrate it into the build process and the source code.

### Adding to the build process

Simply add all of your class files to the `CMakeLists.txt` in the `/SPlisHSPlasH/` directory. We suggest creating new variables for the header and source files and adding these to the `add_library()` as well as to new `source_group()` calls. A possible implementation following our class file conventions would look like the following:

```cmake
set(MYPRESSURESOLVER_HEADER_FILES
    MyPressureSolver/SimulationDataMyPressureSolver.h
    MyPressureSolver/TimeStepMyPressureSolver.h
)

set(MYPRESSURESOLVER_SOURCE_FILES
    MyPressureSolver/SimulationDataMyPressureSolver.cpp
    MyPressureSolver/TimeStepMyPressureSolver.cpp
)

add_library(SPlisHPlasH

    [...]

    ${MYPRESSURESOLVER_HEADER_FILES}
    ${MYPRESSURESOLVER_SOURCE_FILES}
)

source_group("Header Files\\MyPressureSolver" FILES ${MYPRESSURESOLVER_HEADER_FILES})
source_group("Source Files\\MyPressureSolver" FILES ${MYPRESSURESOLVER_SOURCE_FILES})
```

### Integration in the source code

Any timestep method and thus any pressure solver is registered in the `Simulation.h` and `Simulation.cpp` files, which can be found in the `/SPlisHSPlasH/` directory. Adding a new method comprises of the following steps:

* Adding a new enum in `SimulationMethods`
* Creating a new static variable `static int ENUM_SIMULATION_MYPRESSURESOLVER` for the GenericParameter system and initializing it in `Simulation.cpp`
* Including `SPlisHSPlasH/MyPressureSolver/TimeStepMyPressureSolver.h` in `Simulation.cpp`
* Adding a new enum value for `SIMULATION_METHOD` inside `Simulation::initParameters()` using the following line:

```cpp
enumParam->addEnumValue("MyPressureSolverName", ENUM_SIMULATION_MYPRESSURESOLVER);
```

* Adding the pressure solver to `Simulation::setSimulationMethod(...)`, thus making it available for the simulation using the following:

```cpp
else if (method == SimulationMethods::MyPressureSolver)
{
    m_timeStep = new TimeStepMyPressureSolver();
    m_timeStep->init();
    setValue(Simulation::KERNEL_METHOD, <desired standard SPH kernel>);
    setValue(Simulation::GRAD_KERNEL_METHOD, <desired standard SPH gradient kernel>);
}
```

After these additions and building SPlisHSPlasH, our new pressure solver is available inside the simulation.