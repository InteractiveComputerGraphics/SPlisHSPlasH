#include "Simulation.h"
#include "TimeManager.h"
#include "Utilities/Timing.h"
#include "TimeStep.h"
#include "EmitterSystem.h"
#include "SPlisHSPlasH/WCSPH/TimeStepWCSPH.h"
#include "SPlisHSPlasH/PCISPH/TimeStepPCISPH.h"
#include "SPlisHSPlasH/PBF/TimeStepPBF.h"
#include "SPlisHSPlasH/IISPH/TimeStepIISPH.h"
#include "SPlisHSPlasH/DFSPH/TimeStepDFSPH.h"
#include "SPlisHSPlasH/PF/TimeStepPF.h"
#include "BoundaryModel_Akinci2012.h"
#include "BoundaryModel_Bender2019.h"
#include "BoundaryModel_Koschier2017.h"



using namespace SPH;
using namespace std;
using namespace GenParam;

Simulation* Simulation::current = nullptr;
int Simulation::SIM_2D = -1;
int Simulation::PARTICLE_RADIUS = -1;
int Simulation::GRAVITATION = -1;
int Simulation::CFL_METHOD = -1;
int Simulation::CFL_FACTOR = -1;
int Simulation::CFL_MAX_TIMESTEPSIZE = -1;
int Simulation::ENABLE_Z_SORT = -1;
int Simulation::KERNEL_METHOD = -1;
int Simulation::GRAD_KERNEL_METHOD = -1;
int Simulation::ENUM_KERNEL_CUBIC = -1;
int Simulation::ENUM_KERNEL_WENDLANDQUINTICC2 = -1;
int Simulation::ENUM_KERNEL_POLY6 = -1;
int Simulation::ENUM_KERNEL_SPIKY = -1;
int Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC = -1;
int Simulation::ENUM_KERNEL_CUBIC_2D = -1;
int Simulation::ENUM_KERNEL_WENDLANDQUINTICC2_2D = -1;
int Simulation::ENUM_GRADKERNEL_CUBIC = -1;
int Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2 = -1;
int Simulation::ENUM_GRADKERNEL_POLY6 = -1;
int Simulation::ENUM_GRADKERNEL_SPIKY = -1;
int Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC = -1;
int Simulation::ENUM_GRADKERNEL_CUBIC_2D = -1;
int Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D = -1;
int Simulation::SIMULATION_METHOD = -1;
int Simulation::ENUM_CFL_NONE = -1;
int Simulation::ENUM_CFL_STANDARD = -1;
int Simulation::ENUM_CFL_ITER = -1;
int Simulation::ENUM_SIMULATION_WCSPH = -1;
int Simulation::ENUM_SIMULATION_PCISPH = -1;
int Simulation::ENUM_SIMULATION_PBF = -1;
int Simulation::ENUM_SIMULATION_IISPH = -1;
int Simulation::ENUM_SIMULATION_DFSPH = -1;
int Simulation::ENUM_SIMULATION_PF = -1;
int Simulation::BOUNDARY_HANDLING_METHOD = -1;
int Simulation::ENUM_AKINCI2012 = -1;
int Simulation::ENUM_KOSCHIER2017 = -1;
int Simulation::ENUM_BENDER2019 = -1;


Simulation::Simulation () 
{
	m_cflMethod = 1;
	m_cflFactor = 0.5;
	m_cflMaxTimeStepSize = 0.005;
	m_gravitation = Vector3r(0.0, -9.81, 0.0);

	m_kernelMethod = -1;
	m_gradKernelMethod = -1;

	m_neighborhoodSearch = nullptr;
	m_timeStep = nullptr;
	m_simulationMethod = SimulationMethods::NumSimulationMethods;
	m_simulationMethodChanged = NULL;

	m_sim2D = false;
	m_enableZSort = true;

	m_animationFieldSystem = new AnimationFieldSystem();
	m_boundaryHandlingMethod = static_cast<int>(BoundaryHandlingMethods::Bender2019);
}

Simulation::~Simulation () 
{
	delete m_animationFieldSystem;
	delete m_timeStep;
	delete m_neighborhoodSearch;
	delete TimeManager::getCurrent();

	for (unsigned int i = 0; i < m_fluidModels.size(); i++)
		delete m_fluidModels[i];
	m_fluidModels.clear();

	for (unsigned int i = 0; i < m_boundaryModels.size(); i++)
		delete m_boundaryModels[i];
	m_boundaryModels.clear();

	current = nullptr;
}

Simulation* Simulation::getCurrent ()
{
	if (current == nullptr)
	{
		current = new Simulation ();
	}
	return current;
}

void Simulation::setCurrent (Simulation* tm)
{
	current = tm;
}

bool Simulation::hasCurrent()
{
	return (current != nullptr);
}

void Simulation::init(const Real particleRadius, const bool sim2D)
{
	m_sim2D = sim2D;
	initParameters();

	// init kernel
	setParticleRadius(particleRadius);

	setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
	setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);

	// Initialize neighborhood search
	if (m_neighborhoodSearch == NULL)
#ifdef GPU_NEIGHBORHOOD_SEARCH
		m_neighborhoodSearch = new NeighborhoodSearch(m_supportRadius);
#else
		m_neighborhoodSearch = new NeighborhoodSearch(m_supportRadius, false);
#endif
	m_neighborhoodSearch->set_radius(m_supportRadius);
}

void Simulation::initParameters()
{
	ParameterObject::initParameters();

	SIM_2D = createBoolParameter("sim2D", "2D Simulation", &m_sim2D);
	setGroup(SIM_2D, "Simulation");
	setDescription(SIM_2D, "2D/3D simulation.");
	getParameter(SIM_2D)->setReadOnly(true);

	ENABLE_Z_SORT = createBoolParameter("enableZSort", "Enable z-sort", &m_enableZSort);
	setGroup(ENABLE_Z_SORT, "Simulation");
	setDescription(ENABLE_Z_SORT, "Enable z-sort to improve cache hits.");

	ParameterBase::GetFunc<Real> getRadiusFct = std::bind(&Simulation::getParticleRadius, this);
	ParameterBase::SetFunc<Real> setRadiusFct = std::bind(&Simulation::setParticleRadius, this, std::placeholders::_1);
	PARTICLE_RADIUS = createNumericParameter("particleRadius", "Particle radius", getRadiusFct, setRadiusFct);
	setGroup(PARTICLE_RADIUS, "Simulation");
	setDescription(PARTICLE_RADIUS, "Radius of the fluid particles.");
	getParameter(PARTICLE_RADIUS)->setReadOnly(true);

 	GRAVITATION = createVectorParameter("gravitation", "Gravitation", 3u, m_gravitation.data());
 	setGroup(GRAVITATION, "Simulation");
 	setDescription(GRAVITATION, "Vector to define the gravitational acceleration.");

	CFL_METHOD = createEnumParameter("cflMethod", "CFL - method", &m_cflMethod);
	setGroup(CFL_METHOD, "CFL");
	setDescription(CFL_METHOD, "CFL method used for adaptive time stepping.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(CFL_METHOD));
	enumParam->addEnumValue("None", ENUM_CFL_NONE);
	enumParam->addEnumValue("CFL", ENUM_CFL_STANDARD);
	enumParam->addEnumValue("CFL - iterations", ENUM_CFL_ITER);

	CFL_FACTOR = createNumericParameter("cflFactor", "CFL - factor", &m_cflFactor);
	setGroup(CFL_FACTOR, "CFL");
	setDescription(CFL_FACTOR, "Factor to scale the CFL time step size.");
	static_cast<RealParameter*>(getParameter(CFL_FACTOR))->setMinValue(1e-6);

	CFL_MAX_TIMESTEPSIZE = createNumericParameter("cflMaxTimeStepSize", "CFL - max. time step size", &m_cflMaxTimeStepSize);
	setGroup(CFL_MAX_TIMESTEPSIZE, "CFL");
	setDescription(CFL_MAX_TIMESTEPSIZE, "Max. time step size.");
	static_cast<RealParameter*>(getParameter(CFL_MAX_TIMESTEPSIZE))->setMinValue(1e-6);

	ParameterBase::GetFunc<int> getKernelFct = std::bind(&Simulation::getKernel, this);
	ParameterBase::SetFunc<int> setKernelFct = std::bind(&Simulation::setKernel, this, std::placeholders::_1);
	KERNEL_METHOD = createEnumParameter("kernel", "Kernel", getKernelFct, setKernelFct);
	setGroup(KERNEL_METHOD, "Kernel");
	setDescription(KERNEL_METHOD, "Kernel function used in the SPH model.");
	enumParam = static_cast<EnumParameter*>(getParameter(KERNEL_METHOD));
	if (!m_sim2D)
	{
		enumParam->addEnumValue("Cubic spline", ENUM_KERNEL_CUBIC);
		enumParam->addEnumValue("Wendland quintic C2", ENUM_KERNEL_WENDLANDQUINTICC2);
		enumParam->addEnumValue("Poly6", ENUM_KERNEL_POLY6);
		enumParam->addEnumValue("Spiky", ENUM_KERNEL_SPIKY);
		enumParam->addEnumValue("Precomputed cubic spline", ENUM_KERNEL_PRECOMPUTED_CUBIC);
	}
	else
	{
		enumParam->addEnumValue("Cubic spline 2D", ENUM_KERNEL_CUBIC_2D);
		enumParam->addEnumValue("Wendland quintic C2 2D", ENUM_KERNEL_WENDLANDQUINTICC2_2D);
	}

	ParameterBase::GetFunc<int> getGradKernelFct = std::bind(&Simulation::getGradKernel, this);
	ParameterBase::SetFunc<int> setGradKernelFct = std::bind(&Simulation::setGradKernel, this, std::placeholders::_1);
	GRAD_KERNEL_METHOD = createEnumParameter("gradKernel", "Gradient of kernel", getGradKernelFct, setGradKernelFct);
	setGroup(GRAD_KERNEL_METHOD, "Kernel");
	setDescription(GRAD_KERNEL_METHOD, "Gradient of the kernel function used in the SPH model.");
	enumParam = static_cast<EnumParameter*>(getParameter(GRAD_KERNEL_METHOD));
	if (!m_sim2D)
	{
		enumParam->addEnumValue("Cubic spline", ENUM_GRADKERNEL_CUBIC);
		enumParam->addEnumValue("Wendland quintic C2", ENUM_GRADKERNEL_WENDLANDQUINTICC2);
		enumParam->addEnumValue("Poly6", ENUM_GRADKERNEL_POLY6);
		enumParam->addEnumValue("Spiky", ENUM_GRADKERNEL_SPIKY);
		enumParam->addEnumValue("Precomputed cubic spline", ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else
	{
		enumParam->addEnumValue("Cubic spline 2D", ENUM_GRADKERNEL_CUBIC_2D);
		enumParam->addEnumValue("Wendland quintic C2 2D", ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D);
	}

	ParameterBase::GetFunc<int> getSimulationFct = std::bind(&Simulation::getSimulationMethod, this);
	ParameterBase::SetFunc<int> setSimulationFct = std::bind(&Simulation::setSimulationMethod, this, std::placeholders::_1);
	SIMULATION_METHOD = createEnumParameter("simulationMethod", "Simulation method", getSimulationFct, setSimulationFct);
	setGroup(SIMULATION_METHOD, "Simulation");
	setDescription(SIMULATION_METHOD, "Simulation method.");
	enumParam = static_cast<EnumParameter*>(getParameter(SIMULATION_METHOD));
	enumParam->addEnumValue("WCSPH", ENUM_SIMULATION_WCSPH);
	enumParam->addEnumValue("PCISPH", ENUM_SIMULATION_PCISPH);
	enumParam->addEnumValue("PBF", ENUM_SIMULATION_PBF);
	enumParam->addEnumValue("IISPH", ENUM_SIMULATION_IISPH);
	enumParam->addEnumValue("DFSPH", ENUM_SIMULATION_DFSPH);
	enumParam->addEnumValue("Projective Fluids", ENUM_SIMULATION_PF);

	BOUNDARY_HANDLING_METHOD = createEnumParameter("boundaryHandlingMethod", "Boundary handling method", &m_boundaryHandlingMethod);
	setGroup(BOUNDARY_HANDLING_METHOD, "Simulation");
	setDescription(BOUNDARY_HANDLING_METHOD, "Boundary handling method.");
	enumParam = static_cast<EnumParameter*>(getParameter(BOUNDARY_HANDLING_METHOD));
	enumParam->addEnumValue("Akinci et al. 2012", ENUM_AKINCI2012);
	enumParam->addEnumValue("Koschier and Bender 2017", ENUM_KOSCHIER2017);
	enumParam->addEnumValue("Bender et al. 2019", ENUM_BENDER2019);
	enumParam->setReadOnly(true);
}


void Simulation::setParticleRadius(Real val)
{
	m_particleRadius = val;
	m_supportRadius = static_cast<Real>(4.0)*m_particleRadius;

	// init kernel
	Poly6Kernel::setRadius(m_supportRadius);
	SpikyKernel::setRadius(m_supportRadius);
	CubicKernel::setRadius(m_supportRadius);
	WendlandQuinticC2Kernel::setRadius(m_supportRadius);
	PrecomputedCubicKernel::setRadius(m_supportRadius);
	CohesionKernel::setRadius(m_supportRadius);
	AdhesionKernel::setRadius(m_supportRadius);
	CubicKernel2D::setRadius(m_supportRadius);
	WendlandQuinticC2Kernel2D::setRadius(m_supportRadius);
}

void Simulation::setGradKernel(int val)
{
	m_gradKernelMethod = val;

	if (!m_sim2D)
	{
		if ((m_gradKernelMethod < 0) || (m_gradKernelMethod > 4))
			m_gradKernelMethod = 0;

		if (m_gradKernelMethod == 0)
			m_gradKernelFct = CubicKernel::gradW;
		else if (m_gradKernelMethod == 1)
			m_gradKernelFct = WendlandQuinticC2Kernel::gradW;
		else if (m_gradKernelMethod == 2)
			m_gradKernelFct = Poly6Kernel::gradW;
		else if (m_gradKernelMethod == 3)
			m_gradKernelFct = SpikyKernel::gradW;
		else if (m_gradKernelMethod == 4)
			m_gradKernelFct = Simulation::PrecomputedCubicKernel::gradW;
	}
	else
	{
		if ((m_gradKernelMethod < 0) || (m_gradKernelMethod > 1))
			m_gradKernelMethod = 0;

		if (m_gradKernelMethod == 0)
			m_gradKernelFct = CubicKernel2D::gradW;
		else if (m_gradKernelMethod == 1)
			m_gradKernelFct = WendlandQuinticC2Kernel2D::gradW;
	}
}

void Simulation::setKernel(int val)
{
	if (val == m_kernelMethod)
		return;

	m_kernelMethod = val;
	if (!m_sim2D)
	{
		if ((m_kernelMethod < 0) || (m_kernelMethod > 4))
			m_kernelMethod = 0;

		if (m_kernelMethod == 0)
		{
			m_W_zero = CubicKernel::W_zero();
			m_kernelFct = CubicKernel::W;
		}
		else if (m_kernelMethod == 1)
		{
			m_W_zero = WendlandQuinticC2Kernel::W_zero();
			m_kernelFct = WendlandQuinticC2Kernel::W;
		}
		else if (m_kernelMethod == 2)
		{
			m_W_zero = Poly6Kernel::W_zero();
			m_kernelFct = Poly6Kernel::W;
		}
		else if (m_kernelMethod == 3)
		{
			m_W_zero = SpikyKernel::W_zero();
			m_kernelFct = SpikyKernel::W;
		}
		else if (m_kernelMethod == 4)
		{
			m_W_zero = Simulation::PrecomputedCubicKernel::W_zero();
			m_kernelFct = Simulation::PrecomputedCubicKernel::W;
		}
	}
	else
	{
		if ((m_kernelMethod < 0) || (m_kernelMethod > 1))
			m_kernelMethod = 0;

		if (m_kernelMethod == 0)
		{
			m_W_zero = CubicKernel2D::W_zero();
			m_kernelFct = CubicKernel2D::W;
		}
		else if (m_kernelMethod == 1)
		{
			m_W_zero = WendlandQuinticC2Kernel2D::W_zero();
			m_kernelFct = WendlandQuinticC2Kernel2D::W;
		}
	}
	if (getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		updateBoundaryVolume();
}

void Simulation::updateTimeStepSize()
{
	if (m_cflMethod == 1)
		updateTimeStepSizeCFL(0.0001);
	else if (m_cflMethod == 2)
	{
		Real h = TimeManager::getCurrent()->getTimeStepSize();
		updateTimeStepSizeCFL(0.0001);
		const unsigned int iterations = m_timeStep->getValue<unsigned int>(TimeStep::SOLVER_ITERATIONS);
		if (iterations > 10)
			h *= 0.9;
		else if (iterations < 5)
			h *= 1.1;
		h = min(h, TimeManager::getCurrent()->getTimeStepSize());
		TimeManager::getCurrent()->setTimeStepSize(h);
	}
}

void Simulation::updateTimeStepSizeCFL(const Real minTimeStepSize)
{
	const Real radius = m_particleRadius;
	Real h = TimeManager::getCurrent()->getTimeStepSize();
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	// Approximate max. position change due to current velocities
	Real maxVel = 0.1;
	const Real diameter = static_cast<Real>(2.0)*radius;

	// fluid particles
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < numberOfFluidModels(); fluidModelIndex++)
	{
		FluidModel *fm = getFluidModel(fluidModelIndex);
		const unsigned int numParticles = fm->numActiveParticles();
		for (unsigned int i = 0; i < numParticles; i++)
		{
			const Vector3r &vel = fm->getVelocity(i);
			const Vector3r &accel = fm->getAcceleration(i);
			const Real velMag = (vel + accel*h).squaredNorm();
			if (velMag > maxVel)
				maxVel = velMag;
		}
	}

	// boundary particles
	if (getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		{
			BoundaryModel_Akinci2012 *bm = static_cast<BoundaryModel_Akinci2012*>(getBoundaryModel(i));
			if (bm->getRigidBodyObject()->isDynamic())
			{
				for (unsigned int j = 0; j < bm->numberOfParticles(); j++)
				{
					const Vector3r &vel = bm->getVelocity(j);
					const Real velMag = vel.squaredNorm();
					if (velMag > maxVel)
						maxVel = velMag;
				}
			}
		}
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		for (unsigned int boundaryModelIndex = 0; boundaryModelIndex < numberOfBoundaryModels(); boundaryModelIndex++)
		{
			BoundaryModel_Koschier2017 *bm = static_cast<BoundaryModel_Koschier2017*>(getBoundaryModel(boundaryModelIndex));
			if (bm->getRigidBodyObject()->isDynamic())
			{
				maxVel = std::max(maxVel, bm->getMaxVel());
			}
		}
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		for (unsigned int boundaryModelIndex = 0; boundaryModelIndex < numberOfBoundaryModels(); boundaryModelIndex++)
		{
			BoundaryModel_Bender2019 *bm = static_cast<BoundaryModel_Bender2019*>(getBoundaryModel(boundaryModelIndex));
			if (bm->getRigidBodyObject()->isDynamic())
			{
				maxVel = std::max(maxVel, bm->getMaxVel());
			}
		}
	}

	// Approximate max. time step size 		
	h = m_cflFactor * static_cast<Real>(0.4) * (diameter / (sqrt(maxVel)));

	h = min(h, m_cflMaxTimeStepSize);
	h = max(h, minTimeStepSize);

	TimeManager::getCurrent()->setTimeStepSize(h);
}

void Simulation::computeNonPressureForces()
{
	START_TIMING("computeNonPressureForces")
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->computeSurfaceTension();
		fm->computeViscosity();
		fm->computeVorticity();
		fm->computeDragForce();
		fm->computeElasticity();
	}
	STOP_TIMING_AVG
}

void Simulation::reset()
{
	// reset fluid models
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
		getFluidModel(i)->reset();

	// reset boundary models
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		getBoundaryModel(i)->reset();

	if (getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		updateBoundaryVolume();

	if (m_timeStep)
		m_timeStep->reset();

	m_animationFieldSystem->reset();

	performNeighborhoodSearchSort();

	TimeManager::getCurrent()->setTime(0.0);
}

void Simulation::setSimulationMethod(const int val)
{
	SimulationMethods method = static_cast<SimulationMethods>(val);
	if ((method < SimulationMethods::WCSPH) || (method >= SimulationMethods::NumSimulationMethods))
		method = SimulationMethods::DFSPH;

	if (method == m_simulationMethod)
		return;

	delete m_timeStep;
	m_timeStep = nullptr;

	m_simulationMethod = method;

	if (method == SimulationMethods::WCSPH)
	{
		m_timeStep = new TimeStepWCSPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::PCISPH)
	{
		m_timeStep = new TimeStepPCISPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::PBF)
	{
		m_timeStep = new TimeStepPBF();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_POLY6);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_SPIKY);
	}
	else if (method == SimulationMethods::IISPH)
	{
		m_timeStep = new TimeStepIISPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::DFSPH)
	{
		m_timeStep = new TimeStepDFSPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else if (method == SimulationMethods::PF)
	{
		m_timeStep = new TimeStepPF();
		m_timeStep->init();

		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}

	if (m_simulationMethodChanged != nullptr)
		m_simulationMethodChanged();
}



void Simulation::performNeighborhoodSearch()
{
	START_TIMING("neighborhood_search");
	m_neighborhoodSearch->find_neighbors();
	STOP_TIMING_AVG;
}

void Simulation::performNeighborhoodSearchSort()
{
	m_neighborhoodSearch->z_sort();

	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->performNeighborhoodSearchSort();
	}
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
	{
		BoundaryModel *bm = getBoundaryModel(i);
		bm->performNeighborhoodSearchSort();
	}
}

void Simulation::setSimulationMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_simulationMethodChanged = callBackFct;
}

void Simulation::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	model->emittedParticles(startIndex);
	m_timeStep->emittedParticles(model, startIndex);
}

void Simulation::emitParticles()
{
	START_TIMING("emitParticles");
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->getEmitterSystem()->step();
	}
	STOP_TIMING_AVG
}

void Simulation::animateParticles()
{
	START_TIMING("animateParticles");
	m_animationFieldSystem->step();
	STOP_TIMING_AVG
}

void Simulation::addBoundaryModel(BoundaryModel *bm)
{
	m_boundaryModels.push_back(bm);
}

void Simulation::addFluidModel(const std::string &id, const unsigned int nFluidParticles, Vector3r* fluidParticles, Vector3r* fluidVelocities, const unsigned int nMaxEmitterParticles)
{
	FluidModel *fm = new FluidModel();
	fm->initModel(id, nFluidParticles, fluidParticles, fluidVelocities, nMaxEmitterParticles);
	m_fluidModels.push_back(fm);
}


void Simulation::updateBoundaryVolume()
{
	if (m_neighborhoodSearch == nullptr)
		return;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute value psi for boundary particles (boundary handling)
	// (see Akinci et al. "Versatile rigid - fluid coupling for incompressible SPH", Siggraph 2012
	//////////////////////////////////////////////////////////////////////////

	// Search boundary neighborhood

	// Activate only static boundaries
	LOG_INFO << "Initialize boundary volume";
	m_neighborhoodSearch->set_active(false);
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
	{
		if (!getBoundaryModel(i)->getRigidBodyObject()->isDynamic())
			m_neighborhoodSearch->set_active(i + nFluids, true, true);
	}

	//performNeighborhoodSearchSort();
	m_neighborhoodSearch->find_neighbors();

	// Boundary objects
	for (unsigned int body = 0; body < numberOfBoundaryModels(); body++)
	{
		if (!getBoundaryModel(body)->getRigidBodyObject()->isDynamic())
			static_cast<BoundaryModel_Akinci2012*>(getBoundaryModel(body))->computeBoundaryVolume();
	}

	////////////////////////////////////////////////////////////////////////// 
	// Compute boundary psi for all dynamic bodies
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int body = 0; body < numberOfBoundaryModels(); body++)
	{
		// Deactivate all
		m_neighborhoodSearch->set_active(false);

		// Only activate next dynamic body
		if (getBoundaryModel(body)->getRigidBodyObject()->isDynamic())
		{
			m_neighborhoodSearch->set_active(body + nFluids, true, true);
			m_neighborhoodSearch->find_neighbors();
			static_cast<BoundaryModel_Akinci2012*>(getBoundaryModel(body))->computeBoundaryVolume();
		}
	}

	// Activate only fluids 
	m_neighborhoodSearch->set_active(false);
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		for (unsigned int j = 0; j < numberOfFluidModels(); j++)
			m_neighborhoodSearch->set_active(i, j, true);
		for (unsigned int j = numberOfFluidModels(); j < m_neighborhoodSearch->point_sets().size(); j++)
			m_neighborhoodSearch->set_active(i, j, true);
	}
}

void SPH::Simulation::saveState(BinaryFileWriter &binWriter)
{
	binWriter.write(m_W_zero);
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
		getFluidModel(i)->saveState(binWriter);
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		getBoundaryModel(i)->saveState(binWriter);
	m_timeStep->saveState(binWriter);
}

void SPH::Simulation::loadState(BinaryFileReader &binReader)
{
	binReader.read(m_W_zero);
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
		getFluidModel(i)->loadState(binReader);
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		getBoundaryModel(i)->loadState(binReader);
	m_timeStep->loadState(binReader);
}
