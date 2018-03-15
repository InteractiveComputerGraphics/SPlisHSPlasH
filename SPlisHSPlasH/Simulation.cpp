#include "Simulation.h"
#include "SurfaceTension/SurfaceTension_Becker2007.h"
#include "SurfaceTension/SurfaceTension_Akinci2013.h"
#include "SurfaceTension/SurfaceTension_He2014.h"
#include "Viscosity/Viscosity_XSPH.h"
#include "Viscosity/Viscosity_Standard.h"
#include "Viscosity/Viscosity_Bender2017.h"
#include "Viscosity/Viscosity_Peer2015.h"
#include "Viscosity/Viscosity_Peer2016.h"
#include "Viscosity/Viscosity_Takahashi2015.h"
#include "Viscosity/Viscosity_Weiler2018.h"
#include "Vorticity/VorticityConfinement.h"
#include "Drag/DragForce_Gissler2017.h"
#include "Drag/DragForce_Macklin2014.h"
#include "TimeManager.h"
#include "Utilities/Timing.h"
#include "Vorticity/MicropolarModel_Bender2017.h"
#include "TimeStep.h"
#include "SPlisHSPlasH/WCSPH/TimeStepWCSPH.h"
#include "SPlisHSPlasH/PCISPH/TimeStepPCISPH.h"
#include "SPlisHSPlasH/PBF/TimeStepPBF.h"
#include "SPlisHSPlasH/IISPH/TimeStepIISPH.h"
#include "SPlisHSPlasH/DFSPH/TimeStepDFSPH.h"
#include "SPlisHSPlasH/PF/TimeStepPF.h"



using namespace SPH;
using namespace std;
using namespace GenParam;

Simulation* Simulation::current = nullptr;
int Simulation::GRAVITATION = -1;
int Simulation::CFL_METHOD = -1;
int Simulation::CFL_FACTOR = -1;
int Simulation::CFL_MAX_TIMESTEPSIZE = -1;
int Simulation::SIMULATION_METHOD = -1;
int Simulation::DRAG_METHOD = -1;
int Simulation::SURFACE_TENSION_METHOD = -1;
int Simulation::VISCOSITY_METHOD = -1;
int Simulation::VORTICITY_METHOD = -1;
int Simulation::ENUM_CFL_NONE = -1;
int Simulation::ENUM_CFL_STANDARD = -1;
int Simulation::ENUM_CFL_ITER = -1;
int Simulation::ENUM_DRAG_NONE = -1;
int Simulation::ENUM_DRAG_MACKLIN2014 = -1;
int Simulation::ENUM_DRAG_GISSLER2017 = -1;
int Simulation::ENUM_SURFACETENSION_NONE = -1;
int Simulation::ENUM_SURFACETENSION_BECKER2007 = -1;
int Simulation::ENUM_SURFACETENSION_AKINCI2013 = -1;
int Simulation::ENUM_SURFACETENSION_HE2014 = -1;
int Simulation::ENUM_VISCOSITY_NONE = -1;
int Simulation::ENUM_VISCOSITY_STANDARD = -1;
int Simulation::ENUM_VISCOSITY_XSPH = -1;
int Simulation::ENUM_VISCOSITY_BENDER2017 = -1;
int Simulation::ENUM_VISCOSITY_PEER2015 = -1;
int Simulation::ENUM_VISCOSITY_PEER2016 = -1;
int Simulation::ENUM_VISCOSITY_TAKAHASHI2015 = -1;
int Simulation::ENUM_VISCOSITY_WEILER2018 = -1;
int Simulation::ENUM_VORTICITY_NONE = -1;
int Simulation::ENUM_VORTICITY_MICROPOLAR = -1;
int Simulation::ENUM_VORTICITY_VC = -1;
int Simulation::ENUM_SIMULATION_WCSPH = -1;
int Simulation::ENUM_SIMULATION_PCISPH = -1;
int Simulation::ENUM_SIMULATION_PBF = -1;
int Simulation::ENUM_SIMULATION_IISPH = -1;
int Simulation::ENUM_SIMULATION_DFSPH = -1;
int Simulation::ENUM_SIMULATION_PF = -1;


Simulation::Simulation () 
{
	m_cflMethod = 1;
	m_cflFactor = 0.5;
	m_cflMaxTimeStepSize = 0.005;
	m_gravitation = Vector3r(0.0, -9.81, 0.0);

	m_timeStep = nullptr;
	m_simulationMethod = SimulationMethods::WCSPH;
	m_viscosity = nullptr;
	m_viscosityMethod = ViscosityMethods::None;
	m_surfaceTension = nullptr;
	m_surfaceTensionMethod = SurfaceTensionMethods::None;
	m_vorticityMethod = VorticityMethods::None;
	m_vorticity = nullptr;
	m_dragMethod = DragMethods::None;
	m_drag = nullptr;
	m_dragMethodChanged = nullptr;
	m_surfaceTensionMethodChanged = nullptr;
	m_viscosityMethodChanged = nullptr;
	m_vorticityMethodChanged = nullptr;
	m_simulationMethodChanged = NULL;

}

Simulation::~Simulation () 
{
	delete m_surfaceTension;
	delete m_drag;
	delete m_vorticity;
	delete m_viscosity;
	delete m_model;
	delete m_timeStep;
	delete TimeManager::getCurrent();

	current = nullptr;
}

Simulation* Simulation::getCurrent ()
{
	if (current == nullptr)
	{
		current = new Simulation ();
		current->init();
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

void Simulation::init()
{
	initParameters();

	m_model = new FluidModel();
	m_model->init();
	setSimulationMethod(static_cast<int>(SimulationMethods::DFSPH));
	setViscosityMethod(static_cast<int>(ViscosityMethods::Standard));
}

void Simulation::initParameters()
{
	ParameterObject::initParameters();

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

	ParameterBase::GetFunc<int> getDragFct = std::bind(&Simulation::getDragMethod, this);
	ParameterBase::SetFunc<int> setDragFct = std::bind(&Simulation::setDragMethod, this, std::placeholders::_1);
	DRAG_METHOD = createEnumParameter("dragMethod", "Drag method", getDragFct, setDragFct);
	setGroup(DRAG_METHOD, "Drag force");
	setDescription(DRAG_METHOD, "Method to compute drag forces.");
	enumParam = static_cast<EnumParameter*>(getParameter(DRAG_METHOD));
	enumParam->addEnumValue("None", ENUM_DRAG_NONE);
	enumParam->addEnumValue("Macklin et al. 2014", ENUM_DRAG_MACKLIN2014);
	enumParam->addEnumValue("Gissler et al. 2017", ENUM_DRAG_GISSLER2017);


	ParameterBase::GetFunc<int> getSurfaceTensionFct = std::bind(&Simulation::getSurfaceTensionMethod, this);
	ParameterBase::SetFunc<int> setSurfaceTensionFct = std::bind(&Simulation::setSurfaceTensionMethod, this, std::placeholders::_1);
	SURFACE_TENSION_METHOD = createEnumParameter("surfaceTensionMethod", "Surface tension", getSurfaceTensionFct, setSurfaceTensionFct);
	setGroup(SURFACE_TENSION_METHOD, "Surface tension");
	setDescription(SURFACE_TENSION_METHOD, "Method to compute surface tension forces.");
	enumParam = static_cast<EnumParameter*>(getParameter(SURFACE_TENSION_METHOD));
	enumParam->addEnumValue("None", ENUM_SURFACETENSION_NONE);
	enumParam->addEnumValue("Becker & Teschner 2007", ENUM_SURFACETENSION_BECKER2007);
	enumParam->addEnumValue("Akinci et al. 2013", ENUM_SURFACETENSION_AKINCI2013);
	enumParam->addEnumValue("He et al. 2014", ENUM_SURFACETENSION_HE2014);


	ParameterBase::GetFunc<int> getViscosityFct = std::bind(&Simulation::getViscosityMethod, this);
	ParameterBase::SetFunc<int> setViscosityFct = std::bind(&Simulation::setViscosityMethod, this, std::placeholders::_1);
	VISCOSITY_METHOD = createEnumParameter("viscosityMethod", "Viscosity", getViscosityFct, setViscosityFct);
	setGroup(VISCOSITY_METHOD, "Viscosity");
	setDescription(VISCOSITY_METHOD, "Method to compute viscosity forces.");
	enumParam = static_cast<EnumParameter*>(getParameter(VISCOSITY_METHOD));
	enumParam->addEnumValue("None", ENUM_VISCOSITY_NONE);
	enumParam->addEnumValue("Standard", ENUM_VISCOSITY_STANDARD);
	enumParam->addEnumValue("XSPH", ENUM_VISCOSITY_XSPH);
	enumParam->addEnumValue("Bender and Koschier 2017", ENUM_VISCOSITY_BENDER2017);
	enumParam->addEnumValue("Peer et al. 2015", ENUM_VISCOSITY_PEER2015);
	enumParam->addEnumValue("Peer et al. 2016", ENUM_VISCOSITY_PEER2016);
	enumParam->addEnumValue("Takahashi et al. 2015 (improved)", ENUM_VISCOSITY_TAKAHASHI2015);
	enumParam->addEnumValue("Weiler et al. 2018", ENUM_VISCOSITY_WEILER2018);


	ParameterBase::GetFunc<int> getVorticityFct = std::bind(&Simulation::getVorticityMethod, this);
	ParameterBase::SetFunc<int> setVorticityFct = std::bind(&Simulation::setVorticityMethod, this, std::placeholders::_1);
	VORTICITY_METHOD = createEnumParameter("vorticityMethod", "Vorticity", getVorticityFct, setVorticityFct);
	setGroup(VORTICITY_METHOD, "Vorticity");
	setDescription(VORTICITY_METHOD, "Method to compute vorticity forces.");
	enumParam = static_cast<EnumParameter*>(getParameter(VORTICITY_METHOD));
	enumParam->addEnumValue("None", ENUM_VORTICITY_NONE);
	enumParam->addEnumValue("Micropolar model", ENUM_VORTICITY_MICROPOLAR);
	enumParam->addEnumValue("Vorticity confinement", ENUM_VORTICITY_VC);
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
	const Real radius = m_model->getValue<Real>(FluidModel::PARTICLE_RADIUS);
	Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Approximate max. position change due to current velocities
	Real maxVel = 0.1;
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real diameter = 2.0*radius;
	for (unsigned int i = 0; i < numParticles; i++)
	{
		const Vector3r &vel = m_model->getVelocity(0, i);
		const Vector3r &accel = m_model->getAcceleration(i);
		const Real velMag = (vel + accel*h).squaredNorm();
		if (velMag > maxVel)
			maxVel = velMag;
	}

	// boundary particles
	for (unsigned int i = 0; i < m_model->numberOfRigidBodyParticleObjects(); i++)
	{
		FluidModel::RigidBodyParticleObject *rbpo = m_model->getRigidBodyParticleObject(i);
		if (rbpo->m_rigidBody->isDynamic())
		{
			for (unsigned int j = 0; j < rbpo->numberOfParticles(); j++)
			{
				const Vector3r &vel = rbpo->m_v[j];
				const Real velMag = vel.squaredNorm();
				if (velMag > maxVel)
					maxVel = velMag;
			}
		}
	}

	// Approximate max. time step size 		
	h = m_cflFactor * .4 * (diameter / (sqrt(maxVel)));

	h = min(h, m_cflMaxTimeStepSize);
	h = max(h, minTimeStepSize);

	TimeManager::getCurrent()->setTimeStepSize(h);
}

void Simulation::computeNonPressureForces()
{
	START_TIMING("computeNonPressureForces")
	computeSurfaceTension();
	computeViscosity();
	computeVorticity();
	computeDragForce();
	STOP_TIMING_AVG
}

void Simulation::computeSurfaceTension()
{
	if (m_surfaceTension)
		m_surfaceTension->step();
}

void Simulation::computeViscosity()
{
	if (m_viscosity)
		m_viscosity->step();
}

void Simulation::computeVorticity()
{
	if (m_vorticity)
		m_vorticity->step();
}

void Simulation::computeDragForce()
{
	if (m_drag)
		m_drag->step();
}

void Simulation::reset()
{
	m_model->reset();
	if (m_timeStep)
		m_timeStep->reset();
	if (m_surfaceTension)
		m_surfaceTension->reset();
	if (m_viscosity)
		m_viscosity->reset();
	if (m_vorticity)
		m_vorticity->reset();
	if (m_drag)
		m_drag->reset();

	TimeManager::getCurrent()->setTime(0.0);
	TimeManager::getCurrent()->setTimeStepSize(0.001);
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
		setValue(Simulation::CFL_METHOD, Simulation::ENUM_CFL_NONE);
		m_model->setValue(FluidModel::KERNEL_METHOD, FluidModel::ENUM_KERNEL_CUBIC);
		m_model->setValue(FluidModel::GRAD_KERNEL_METHOD, FluidModel::ENUM_GRADKERNEL_CUBIC);
		TimeManager::getCurrent()->setTimeStepSize(0.001);
	}
	else if (method == SimulationMethods::PCISPH)
	{
		m_timeStep = new TimeStepPCISPH();
		m_timeStep->init();
		m_model->setValue(FluidModel::KERNEL_METHOD, FluidModel::ENUM_KERNEL_CUBIC);
		m_model->setValue(FluidModel::GRAD_KERNEL_METHOD, FluidModel::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::PBF)
	{
		m_timeStep = new TimeStepPBF();
		m_timeStep->init();
		m_model->setValue(FluidModel::KERNEL_METHOD, FluidModel::ENUM_KERNEL_POLY6);
		m_model->setValue(FluidModel::GRAD_KERNEL_METHOD, FluidModel::ENUM_GRADKERNEL_SPIKY);
	}
	else if (method == SimulationMethods::IISPH)
	{
		m_timeStep = new TimeStepIISPH();
		m_timeStep->init();
		m_model->setValue(FluidModel::KERNEL_METHOD, FluidModel::ENUM_KERNEL_CUBIC);
		m_model->setValue(FluidModel::GRAD_KERNEL_METHOD, FluidModel::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::DFSPH)
	{
		m_timeStep = new TimeStepDFSPH();
		m_timeStep->init();
		m_model->setValue(FluidModel::KERNEL_METHOD, FluidModel::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		m_model->setValue(FluidModel::GRAD_KERNEL_METHOD, FluidModel::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else if (method == SimulationMethods::PF)
	{
		m_timeStep = new TimeStepPF();
		m_timeStep->init();
		m_model->setValue(FluidModel::KERNEL_METHOD, FluidModel::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		m_model->setValue(FluidModel::GRAD_KERNEL_METHOD, FluidModel::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}

	if (m_simulationMethodChanged != nullptr)
		m_simulationMethodChanged();
}

void Simulation::setSurfaceTensionMethod(const int val)
{
	SurfaceTensionMethods stm = static_cast<SurfaceTensionMethods>(val);
	if ((stm < SurfaceTensionMethods::None) || (stm >= SurfaceTensionMethods::NumSurfaceTensionMethods))
		stm = SurfaceTensionMethods::None;
	if (stm == m_surfaceTensionMethod)
		return;

	delete m_surfaceTension;
	m_surfaceTension = nullptr;

	m_surfaceTensionMethod = stm;
	if (m_surfaceTensionMethod == SurfaceTensionMethods::Becker2007)
		m_surfaceTension = new SurfaceTension_Becker2007(m_model);
	else if (m_surfaceTensionMethod == SurfaceTensionMethods::Akinci2013)
		m_surfaceTension = new SurfaceTension_Akinci2013(m_model);
	else if (m_surfaceTensionMethod == SurfaceTensionMethods::He2014)
		m_surfaceTension = new SurfaceTension_He2014(m_model);

	if (m_surfaceTension != nullptr)
		m_surfaceTension->init();

	if (m_surfaceTensionMethodChanged != nullptr)
		m_surfaceTensionMethodChanged();
}

void Simulation::setViscosityMethod(const int val)
{
	ViscosityMethods vm = static_cast<ViscosityMethods>(val);
	if ((vm < ViscosityMethods::None) || (vm >= ViscosityMethods::NumViscosityMethods))
		vm = ViscosityMethods::XSPH;

	if (vm == m_viscosityMethod)
		return;

	delete m_viscosity;
	m_viscosity = nullptr;

	m_viscosityMethod = vm;

	if (m_viscosityMethod == ViscosityMethods::Standard)
		m_viscosity = new Viscosity_Standard(m_model);
	else if (m_viscosityMethod == ViscosityMethods::XSPH)
		m_viscosity = new Viscosity_XSPH(m_model);
	else if (m_viscosityMethod == ViscosityMethods::Bender2017)
		m_viscosity = new Viscosity_Bender2017(m_model);
	else if (m_viscosityMethod == ViscosityMethods::Peer2015)
		m_viscosity = new Viscosity_Peer2015(m_model);
	else if (m_viscosityMethod == ViscosityMethods::Peer2016)
		m_viscosity = new Viscosity_Peer2016(m_model);
	else if (m_viscosityMethod == ViscosityMethods::Takahashi2015)
		m_viscosity = new Viscosity_Takahashi2015(m_model);
	else if (m_viscosityMethod == ViscosityMethods::Weiler2018)
		m_viscosity = new Viscosity_Weiler2018(m_model);

	if (m_viscosity != nullptr)
		m_viscosity->init();

	if (m_viscosityMethodChanged != nullptr)
		m_viscosityMethodChanged();
}


void Simulation::setVorticityMethod(const int val)
{
	VorticityMethods vm = static_cast<VorticityMethods>(val);
	if ((vm < VorticityMethods::None) || (vm >= VorticityMethods::NumVorticityMethods))
		vm = VorticityMethods::None;

	if (vm == m_vorticityMethod)
		return;

	delete m_vorticity;
	m_vorticity = nullptr;

	m_vorticityMethod = vm;

	if (m_vorticityMethod == VorticityMethods::Micropolar)
		m_vorticity = new MicropolarModel_Bender2017(m_model);
	else if (m_vorticityMethod == VorticityMethods::VorticityConfinement)
		m_vorticity = new VorticityConfinement(m_model);

	if (m_vorticity != nullptr)
		m_vorticity->init();

	if (m_vorticityMethodChanged != nullptr)
		m_vorticityMethodChanged();
}

void Simulation::setDragMethod(const int val)
{
	DragMethods dm = static_cast<DragMethods>(val);
	if ((dm < DragMethods::None) || (dm >= DragMethods::NumDragMethods))
		dm = DragMethods::None;

	if (dm == m_dragMethod)
		return;

	delete m_drag;
	m_drag = nullptr;

	m_dragMethod = dm;

	if (m_dragMethod == DragMethods::Gissler2017)
		m_drag = new DragForce_Gissler2017(m_model);
	else if (m_dragMethod == DragMethods::Macklin2014)
		m_drag = new DragForce_Macklin2014(m_model);

	if (m_drag != nullptr)
		m_drag->init();

	if (m_dragMethodChanged != nullptr)
		m_dragMethodChanged();
}

void Simulation::performNeighborhoodSearch()
{
	START_TIMING("neighborhood_search");
	m_model->getNeighborhoodSearch()->find_neighbors();
	STOP_TIMING_AVG;
}

void Simulation::performNeighborhoodSearchSort()
{
	if (m_viscosity)
		m_viscosity->performNeighborhoodSearchSort();
	if (m_surfaceTension)
		m_surfaceTension->performNeighborhoodSearchSort();
	if (m_vorticity)
		m_vorticity->performNeighborhoodSearchSort();
	if (m_drag)
		m_drag->performNeighborhoodSearchSort();
}

void Simulation::setDragMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_dragMethodChanged = callBackFct;
}

void Simulation::setSurfaceMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_surfaceTensionMethodChanged = callBackFct;
}

void Simulation::setViscosityMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_viscosityMethodChanged = callBackFct;
}

void Simulation::setVorticityMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_vorticityMethodChanged = callBackFct;
}

void Simulation::setSimulationMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_simulationMethodChanged = callBackFct;
}

void Simulation::emittedParticles(const unsigned int startIndex)
{
	if (m_viscosity)
		m_viscosity->emittedParticles(startIndex);
	if (m_surfaceTension)
		m_surfaceTension->emittedParticles(startIndex);
	if (m_vorticity)
		m_vorticity->emittedParticles(startIndex);
	if (m_drag)
		m_drag->emittedParticles(startIndex);
}

void Simulation::emitParticles()
{
	getModel()->getEmitterSystem()->step();
}