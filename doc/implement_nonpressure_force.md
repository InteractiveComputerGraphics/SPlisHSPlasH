# Implementing a new non-pressure force method

Non-pressure forces (e.g. viscosity, vorticity, surface tension or drag forces) are all implemented in the same way in SPlisHSPlasH. In the following we explain the implementation of such a method using as example a new viscosity method. 

SPlisHSPlasH organizes the viscosities in `/SPlisHSPlasH/Viscosity/` and thus any changes or additions are intended to take place in this directory. The user can add new viscosity methods by creating new or copying and modifying existing viscosity class files and registering these inside the build system and the source code.

## Creating a new class

If you want to create a new viscosity class from scratch, you should consider reading the doxygen documentation on the `ViscosityBase` class and several of its derived classes. In short, every viscosity method inherits from the base class `ViscosityBase`, which itself inherits from `NonPressureForceBase`. A minimal working derived class would look like this:

**MyViscosity.h**

```cpp
#ifndef __MyViscosity_h__
#define __MyViscosity_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"

namespace SPH 
{
	class MyViscosity : public ViscosityBase
	{
	protected:
		virtual void initParameters();

	public:
		MyViscosity(FluidModel *model);
		virtual ~MyViscosity(void);
		
		static NonPressureForceBase* creator(FluidModel* model) { return new MyViscosity(model); }

		virtual void step();
		virtual void reset();
	};
}

#endif
```

**MyViscosity.cpp**
```cpp
#include "MyViscosity.h"

MyViscosity::MyViscosity(FluidModel *model) :
	ViscosityBase(model)
{
	[...]
}

MyViscosity::~MyViscosity(void)
{
	[...]
}

void MyViscosity::initParameters()
{
	ViscosityBase::initParameters();

	[...]
}

void MyViscosity::step()
{
	[...]
}

void MyViscosity::reset()
{
	[...]
}
```

including the following:

* a constructor with `FluidModel*` as the sole parameter `MyViscosity(FluidModel *model)`
* a `initParameters()` method calling the base class method for parameter setup
* a step function `void step()` called in each timestep for the associated fluid
* a reset function `void reset()` called on every reset of the simulation

### Customizing your class

#### Neighborhood search sort

The user is also free to add and save additional per particle data inside the viscosity method, but has to ensure that these are also included in the *neighborhood search sort*. Sorting is required if the data is used over multiple simulations steps. The neighborhood search performs a z-sort every n steps to improve the number of cache hits. Since all particles are resorted, also their data must be resorted. For this, the user has to override the `performNeighborhoodSearchSort()` method. A minimal example would look like the following:

```cpp
void MyViscosity::performNeighborhoodSearchSort()
{
	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_myParticleViscosityData[0]);
}
```

#### Additional particle fields

For visualization and/or debugging purposes, the user may also want to subject the particle data to SPlisHSPlasH's particle informations. To do this, the user has to add the particle data field to the list of fields inside each `FluidModel`. This can be for example done in the constructor by adding the `addField(const FieldDescription &field)` of the corresponding `FluidModel`. The fields can be used to define the color of a particle, they can be exported to bgeo or ParaView and in the simulator the user can output the field data of the selected particles by pressing "i".

For more information, please refer to the doxygen documentation and maybe take a look at the already existing implementations. Adding a field has the following form:

```cpp
model->addField({ "myFieldName", <FieldType>, <lambda expression returning reference to the data field>}, <save state (boolean)>);
```
Here is an example:
```cpp
model->addField({ "myFieldName", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_myFieldValues[i][0]; }, true });
```
The field name is used in the GUI and when exporting the data. The boolean at the end determines if this field should be stored when the simulation state is saved. This should only be done if the value is not recomputed in each simulation step so that the value of the last step is required. 

Also don't forget to remove the field, when the instance of the viscosity method is destroyed:

```cpp
m_model->removeFiledByName("myFieldName");
```

#### Deferred initialization

The user can override the `deferredInit()` method. This function is called after the simulation scene is loaded and all parameters are initialized. While reading a scene file several parameters can change. The `deferredInit()` function should initialize all values which depend on these parameters.

```cpp
void MyViscosity::deferredInit()
{
	initMyViscosity();
}
```

## Registering the viscosity method

To add our new viscosity method, we have to integrate it into the build process and the source code. 

### Adding to the build process

Simply add the class files `MyViscosity.h` and `MyViscosity.cpp` to the `CMakeLists.txt` in the `/SPlisHSPlasH/` directory. This can be done by adding the relative file paths to the respective variables `VISCOSITY_HEADER_FILES` and `VISCOSITY_SOURCE_FILES`:

```cmake
set(VISCOSITY_HEADER_FILES
	[...]
	Viscosity/MyViscosity.h
)

set(VISCOSITY_SOURCE_FILES
	[...]
	Viscosity/MyViscosity.cpp
)
```

### Integration in the source code

Any non-pressure force method is registered in the file `NonPressureForceRegistration.cpp`, which can be found in the `/SPlisHSPlasH/` directory. Adding our new viscosity method is done by adding the following line to the function `void Simulation::registerNonpressureForces()`:

```cpp 
addViscosityMethod("My viscosity method", MyViscosity::creator);
```

and including `Viscosity/MyViscosity.h`.

After these additions and building SPlisHSPlasH, our new viscosity method is available inside the simulation.