# Macros

SPlisHSPlasH defines useful macros to e.g. iterate over all neighboring particles inside the neighborhood of the current one. These can be found in `Simulation.h`. In the following, we want to give a short overview over these macros. For further information, please refer to the api documentation.

## Looping over fluid neighbors

An essential part of SPH computation is to use the properties of neighboring particles to compute the desired value. SPlisHSPlasH provides macros iterating over every fluid neighbor, which can be used like predefined for-loop constructs. These include the following:

### forall_fluid_neigbors

```cpp
#define forall_fluid_neighbors(code) \
	for (unsigned int pid = 0; pid < nFluids; pid++) \
	{ \
		FluidModel *fm_neighbor = sim->getFluidModelFromPointSet(pid); \
		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++) \
		{ \
			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j); \
			const Vector3r &xj = fm_neighbor->getPosition(neighborIndex); \
			code \
		} \
	} 
```
`forall_fluid_neigbors` loops over every fluid particle (in all fluid phases) in the neighborhood region of the current one. Note that this does **not** include boundary particles. The user can use this macro by writing the desired code inside the brackets. For the usage of most of the macros, some additional variables have to be predefined. These include in this case: 

- `Simulation *sim = Simulation::getCurrent()`, the current simulation instance
- `unsigned int nFluids`, the amount of FluidModel instances
- `unsigned int fluidModelIndex`, the index of the FluidModel of the current particle
- `unsigned int i`, the index of the current particle inside the FluidModel with index fluidModelIndex

Further, this macro also defines certain variables, which can be accessed inside the code given to the macro:

- `unsigned int pid`, the index of the FluidModel of the neighboring particle
- `FluidModel *fm_neighbor`, the FluidModel reference of the neighboring particle
- `const unsigned int neighborIndex`, the particle index of the neighboring particle
- `const Vector3r &xj`, the position of the neighboring particle

Henceforth, we denote the required additional variables by **Requires** and the by the macro defined ones by **Defines**.

### forall_fluid_neighbors_in_same_phase

```cpp
#define forall_fluid_neighbors_in_same_phase(code) \
	for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++) \
	{ \
		const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j); \
		const Vector3r &xj = model->getPosition(neighborIndex); \
		code \
	} 
```
`forall_fluid_neighbors_in_same_phase` loops over every fluid particle in the neighborhood region considering only neighbors from the **same** FluidModel as the current one.
- **Requires**:
  - `Simulation *sim = Simulation::getCurrent()`
  - `unsigned int fluidModelIndex`
  - `unsigned int i`
- **Defines**:
  - `const unsigned int neighborIndex`
  - `const Vector3r &xj`

## Looping over boundaries

### forall_boundary_neighbors

```cpp
#define forall_boundary_neighbors(code) \
    for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) \
    { \
	    BoundaryModel_Akinci2012 *bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid)); \
	    for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++) \
	    { \
	    	const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j); \
	    	const Vector3r &xj = bm_neighbor->getPosition(neighborIndex); \
	    	code \
	    } \
    }
```
`forall_boundary_neighbors` loops over all boundary neighbors casting them to the Akinci 2012 boundary model.
- **Requires**:
  - `Simulation *sim = Simulation::getCurrent()`
  - `unsigned int nFluids`
  - `unsigned int fluidModelIndex`
  - `unsigned int i`
- **Defines**:
  - `unsigned int pid`, the index of the FluidModel associated with the BoundaryModel
  - `BoundaryModel_Akinci2012 *bm_neighbor`, the BoundaryModel reference of the neighboring particle
  - `const unsigned int neighborIndex`, the particle index of the neighboring particle
  - `const Vector3r &xj`, the position of the neigboring particle 

### forall_density_maps

```cpp
#define forall_density_maps(code) \
for (unsigned int pid = 0; pid < nBoundaries; pid++) \
{ \
	BoundaryModel_Koschier2017 *bm_neighbor = static_cast<BoundaryModel_Koschier2017*>(sim->getBoundaryModel(pid)); \
	const Real rho = bm_neighbor->getBoundaryDensity(fluidModelIndex, i); \
	if (rho != 0.0) \
	{ \
		const Vector3r &gradRho = bm_neighbor->getBoundaryDensityGradient(fluidModelIndex, i).cast<Real>(); \
		const Vector3r &xj = bm_neighbor->getBoundaryXj(fluidModelIndex, i); \
		code \
	} \
}
```
`forall_density_maps` loops over all boundary neighbors casting them to the Koschier 2017 boundary model.
- **Requires**:
  - `Simulation *sim = Simulation::getCurrent()`
  - `unsigned int nBoundaries`
  - `unsigned int fluidModelIndex`
  - `unsigned int i`
- **Defines**:
  - `unsigned int pid`
  - `BoundaryModel_Koschier2017 *bm_neighbor`
  - `const Real rho`, the boundary density given by the density map
  - `const Vector3r &gradRho`, the boundary density gradient
  - `const Vector3r &xj`

### forall_volume_maps

```cpp
#define forall_volume_maps(code) \
    for (unsigned int pid = 0; pid < nBoundaries; pid++) \
    { \
    	BoundaryModel_Bender2019 *bm_neighbor = static_cast<BoundaryModel_Bender2019*>(sim->getBoundaryModel(pid)); \
    	const Real Vj = bm_neighbor->getBoundaryVolume(fluidModelIndex, i);  \
    	if (Vj > 0.0) \
    	{ \
    		const Vector3r &xj = bm_neighbor->getBoundaryXj(fluidModelIndex, i); \
    		code \
    	} \
    }
```
`forall_volume_maps` loops over all boundary neighbors casting them to the Bender 2019 boundary model.
- **Requires**:
  - `Simulation *sim = Simulation::getCurrent()`
  - `unsigned int nBoundaries`
  - `unsigned int fluidModelIndex`
  - `unsigned int i`
- **Defines**:
  - `unsigned int pid`
  - `BoundaryModel_Koschier2019 *bm_neighbor`
  - `const Real Vj`, the boundary volume given by the volume map
  - `const Vector3r &xj`

## AVX variants

SPlisHSPlasH also defines versions using AVX optimizations for some of the macros. These can be used if the respective CMake option is set in the building process. Note that many of the aforementioned by the macro defined variables are given in AVX compatible data types, if you choose to use the AVX version of these macros.