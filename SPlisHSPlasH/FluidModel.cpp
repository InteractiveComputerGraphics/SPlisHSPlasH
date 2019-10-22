#include "FluidModel.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "NeighborhoodSearch.h"
#include "Simulation.h"
#include "EmitterSystem.h"
#include "Viscosity/ViscosityBase.h"
#include "SurfaceTension/SurfaceTensionBase.h"
#include "Vorticity/VorticityBase.h"
#include "Drag/DragBase.h"
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
#include "Vorticity/MicropolarModel_Bender2017.h"
#include "Drag/DragForce_Gissler2017.h"
#include "Drag/DragForce_Macklin2014.h"
#include "Elasticity/ElasticityBase.h"
#include "Elasticity/Elasticity_Becker2009.h"
#include "Elasticity/Elasticity_Peer2018.h"


using namespace SPH;
using namespace GenParam;

int FluidModel::NUM_PARTICLES = -1;
int FluidModel::NUM_REUSED_PARTICLES = -1;
int FluidModel::DENSITY0 = -1;
int FluidModel::DRAG_METHOD = -1;
int FluidModel::SURFACE_TENSION_METHOD = -1;
int FluidModel::VISCOSITY_METHOD = -1;
int FluidModel::VORTICITY_METHOD = -1;
int FluidModel::ELASTICITY_METHOD = -1;
int FluidModel::ENUM_DRAG_NONE = -1;
int FluidModel::ENUM_DRAG_MACKLIN2014 = -1;
int FluidModel::ENUM_DRAG_GISSLER2017 = -1;
int FluidModel::ENUM_SURFACETENSION_NONE = -1;
int FluidModel::ENUM_SURFACETENSION_BECKER2007 = -1;
int FluidModel::ENUM_SURFACETENSION_AKINCI2013 = -1;
int FluidModel::ENUM_SURFACETENSION_HE2014 = -1;
int FluidModel::ENUM_VISCOSITY_NONE = -1;
int FluidModel::ENUM_VISCOSITY_STANDARD = -1;
int FluidModel::ENUM_VISCOSITY_XSPH = -1;
int FluidModel::ENUM_VISCOSITY_BENDER2017 = -1;
int FluidModel::ENUM_VISCOSITY_PEER2015 = -1;
int FluidModel::ENUM_VISCOSITY_PEER2016 = -1;
int FluidModel::ENUM_VISCOSITY_TAKAHASHI2015 = -1;
int FluidModel::ENUM_VISCOSITY_WEILER2018 = -1;
int FluidModel::ENUM_VORTICITY_NONE = -1;
int FluidModel::ENUM_VORTICITY_MICROPOLAR = -1;
int FluidModel::ENUM_VORTICITY_VC = -1;
int FluidModel::ENUM_ELASTICITY_NONE = -1;
int FluidModel::ENUM_ELASTICITY_BECKER2009 = -1;
int FluidModel::ENUM_ELASTICITY_PEER2018 = -1;


FluidModel::FluidModel() :
	m_masses(),
	m_a(),
	m_v0(),
	m_x0(),
	m_x(),
	m_v(),
	m_density(),
	m_particleId(),
	m_particleState()
{		
	m_density0 = 1000.0;
	m_pointSetIndex = 0;

	m_emitterSystem = new EmitterSystem(this);
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
	m_elasticityMethod = ElasticityMethods::None;
	m_elasticity = nullptr;
	m_elasticityMethodChanged = nullptr;

	addField({ "id", FieldType::UInt, [&](const unsigned int i) -> unsigned int* { return &getParticleId(i); }, true });
	addField({ "position", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &getPosition(i)[0]; }, true });
	addField({ "velocity", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &getVelocity(i)[0]; }, true });
	addField({ "density", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &getDensity(i); }, false });
}

FluidModel::~FluidModel(void)
{
	removeFieldByName("position");
	removeFieldByName("velocity");
	removeFieldByName("density");

	delete m_emitterSystem;
	delete m_surfaceTension;
	delete m_drag;
	delete m_vorticity;
	delete m_viscosity;
	delete m_elasticity;

	releaseFluidParticles();
}

void FluidModel::init()
{
	initParameters();

	setViscosityMethod(static_cast<int>(ViscosityMethods::Standard));
}

void FluidModel::initParameters()
{
	std::string groupName = "FluidModel";
	ParameterObject::initParameters();

	ParameterBase::GetFunc<Real> getDensity0Fct = std::bind(&FluidModel::getDensity0, this);
	ParameterBase::SetFunc<Real> setDensity0Fct = std::bind(&FluidModel::setDensity0, this, std::placeholders::_1);
	DENSITY0 = createNumericParameter("density0", "Rest density", getDensity0Fct, setDensity0Fct);
	setGroup(DENSITY0, groupName);
	setDescription(DENSITY0, "Rest density of the fluid.");
	getParameter(DENSITY0)->setReadOnly(true);

	NUM_PARTICLES = createNumericParameter("numParticles", "# active particles", &m_numActiveParticles);
	setGroup(NUM_PARTICLES, groupName);
	setDescription(NUM_PARTICLES, "Number of active fluids particles in the simulation.");
	getParameter(NUM_PARTICLES)->setReadOnly(true);

	NUM_REUSED_PARTICLES = createNumericParameter<unsigned int>("numReusedParticles", "# reused particles", [&]() { return m_emitterSystem->numReusedParticles(); });
	setGroup(NUM_REUSED_PARTICLES, groupName);
	setDescription(NUM_REUSED_PARTICLES, "Number of reused fluid particles in the simulation.");
	getParameter(NUM_REUSED_PARTICLES)->setReadOnly(true);

	ParameterBase::GetFunc<int> getDragFct = std::bind(&FluidModel::getDragMethod, this);
	ParameterBase::SetFunc<int> setDragFct = std::bind(&FluidModel::setDragMethod, this, std::placeholders::_1);
	DRAG_METHOD = createEnumParameter("dragMethod", "Drag method", getDragFct, setDragFct);
	setGroup(DRAG_METHOD, "Drag force");
	setDescription(DRAG_METHOD, "Method to compute drag forces.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(DRAG_METHOD));
	enumParam->addEnumValue("None", ENUM_DRAG_NONE);
	enumParam->addEnumValue("Macklin et al. 2014", ENUM_DRAG_MACKLIN2014);
	enumParam->addEnumValue("Gissler et al. 2017", ENUM_DRAG_GISSLER2017);


	ParameterBase::GetFunc<int> getSurfaceTensionFct = std::bind(&FluidModel::getSurfaceTensionMethod, this);
	ParameterBase::SetFunc<int> setSurfaceTensionFct = std::bind(&FluidModel::setSurfaceTensionMethod, this, std::placeholders::_1);
	SURFACE_TENSION_METHOD = createEnumParameter("surfaceTensionMethod", "Surface tension", getSurfaceTensionFct, setSurfaceTensionFct);
	setGroup(SURFACE_TENSION_METHOD, "Surface tension");
	setDescription(SURFACE_TENSION_METHOD, "Method to compute surface tension forces.");
	enumParam = static_cast<EnumParameter*>(getParameter(SURFACE_TENSION_METHOD));
	enumParam->addEnumValue("None", ENUM_SURFACETENSION_NONE);
	enumParam->addEnumValue("Becker & Teschner 2007", ENUM_SURFACETENSION_BECKER2007);
	enumParam->addEnumValue("Akinci et al. 2013", ENUM_SURFACETENSION_AKINCI2013);
	enumParam->addEnumValue("He et al. 2014", ENUM_SURFACETENSION_HE2014);


	ParameterBase::GetFunc<int> getViscosityFct = std::bind(&FluidModel::getViscosityMethod, this);
	ParameterBase::SetFunc<int> setViscosityFct = std::bind(&FluidModel::setViscosityMethod, this, std::placeholders::_1);
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

	ParameterBase::GetFunc<int> getVorticityFct = std::bind(&FluidModel::getVorticityMethod, this);
	ParameterBase::SetFunc<int> setVorticityFct = std::bind(&FluidModel::setVorticityMethod, this, std::placeholders::_1);
	VORTICITY_METHOD = createEnumParameter("vorticityMethod", "Vorticity", getVorticityFct, setVorticityFct);
	setGroup(VORTICITY_METHOD, "Vorticity");
	setDescription(VORTICITY_METHOD, "Method to compute vorticity forces.");
	enumParam = static_cast<EnumParameter*>(getParameter(VORTICITY_METHOD));
	enumParam->addEnumValue("None", ENUM_VORTICITY_NONE);
	enumParam->addEnumValue("Micropolar model", ENUM_VORTICITY_MICROPOLAR);
	enumParam->addEnumValue("Vorticity confinement", ENUM_VORTICITY_VC);

	ParameterBase::GetFunc<int> getElasticityFct = std::bind(&FluidModel::getElasticityMethod, this);
	ParameterBase::SetFunc<int> setElasticityFct = std::bind(&FluidModel::setElasticityMethod, this, std::placeholders::_1);
	ELASTICITY_METHOD = createEnumParameter("elasticityMethod", "Elasticity method", getElasticityFct, setElasticityFct);
	setGroup(ELASTICITY_METHOD, "Elasticity");
	setDescription(ELASTICITY_METHOD, "Method to compute elastic forces.");
	enumParam = static_cast<EnumParameter*>(getParameter(ELASTICITY_METHOD));
	enumParam->addEnumValue("None", ENUM_ELASTICITY_NONE);
	enumParam->addEnumValue("Becker et al. 2009", ENUM_ELASTICITY_BECKER2009);
	enumParam->addEnumValue("Peer et al. 2018", ENUM_ELASTICITY_PEER2018);
}


void FluidModel::reset()
{
	setNumActiveParticles(m_numActiveParticles0);
	const unsigned int nPoints = numActiveParticles();

	// Fluid
	for (unsigned int i = 0; i < nPoints; i++)
	{
		const Vector3r& x0 = getPosition0(i);
		getPosition(i) = x0;
		getVelocity(i) = getVelocity0(i);
		getAcceleration(i).setZero();
		m_density[i] = 0.0;
		m_particleId[i] = i;
		m_particleState[i] = ParticleState::Active;
	}
	// emitted particles
	for (unsigned int i = nPoints; i < (unsigned int)m_particleId.size(); i++)
	{
		m_particleId[i] = i;
	}

	NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
	if (neighborhoodSearch->point_set(m_pointSetIndex).n_points() != nPoints)
		neighborhoodSearch->resize_point_set(m_pointSetIndex, &getPosition(0)[0], nPoints);

	if (m_surfaceTension)
		m_surfaceTension->reset();
	if (m_viscosity)
		m_viscosity->reset();
	if (m_vorticity)
		m_vorticity->reset();
	if (m_drag)
		m_drag->reset();
	if (m_elasticity)
		m_elasticity->reset();

	m_emitterSystem->reset();
}

void FluidModel::initMasses()
{
	const Real particleRadius = Simulation::getCurrent()->getParticleRadius();
	const int nParticles = (int) numParticles();
	const Real diam = static_cast<Real>(2.0)*particleRadius;
	if (Simulation::getCurrent()->is2DSimulation())
		m_V = static_cast<Real>(0.8) * diam*diam;
	else
		m_V = static_cast<Real>(0.8) * diam*diam*diam;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < nParticles; i++)
		{
			setMass(i, m_V * m_density0);						// each particle represents a cube with a side length of r		
																// mass is slightly reduced to prevent pressure at the beginning of the simulation
		}
	}
}

void FluidModel::resizeFluidParticles(const unsigned int newSize)
{
	m_x0.resize(newSize);
	m_x.resize(newSize);
	m_v.resize(newSize);
	m_v0.resize(newSize);
	m_a.resize(newSize);
	m_masses.resize(newSize);
	m_density.resize(newSize);
	m_particleId.resize(newSize);
	m_particleState.resize(newSize);
}

void FluidModel::releaseFluidParticles()
{
	m_x0.clear();
	m_x.clear();
	m_v.clear();
	m_v0.clear();
	m_a.clear();
	m_masses.clear();
	m_density.clear();
	m_particleId.clear();
	m_particleState.clear();
}

void FluidModel::initModel(const std::string &id, const unsigned int nFluidParticles, Vector3r* fluidParticles, Vector3r* fluidVelocities, const unsigned int nMaxEmitterParticles)
{
	m_id = id;
	init();
	releaseFluidParticles();
	resizeFluidParticles(nFluidParticles + nMaxEmitterParticles);

	// copy fluid positions
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)nFluidParticles; i++)
		{
			getPosition0(i) = fluidParticles[i];
			getPosition(i) = fluidParticles[i];
			getVelocity0(i) = fluidVelocities[i];
			getVelocity(i) = fluidVelocities[i];
			getAcceleration(i).setZero();
			m_density[i] = 0.0;
			m_particleId[i] = i;
			m_particleState[i] = ParticleState::Active;
		}
	}
	// set IDs for emitted particles
	for (unsigned int i = nFluidParticles; i < (nFluidParticles + nMaxEmitterParticles); i++)
	{
		m_particleId[i] = i;
	}

	// initialize masses
	initMasses();

	// Fluids 
	NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
	m_pointSetIndex = neighborhoodSearch->add_point_set(&getPosition(0)[0], nFluidParticles, true, true, true, this);

	m_numActiveParticles0 = nFluidParticles;
	m_numActiveParticles = m_numActiveParticles0;
}


void FluidModel::performNeighborhoodSearchSort()
{
	const unsigned int numPart = numActiveParticles();
	if (numPart == 0)
		return;

	NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();

	auto const& d = neighborhoodSearch->point_set(m_pointSetIndex);
	d.sort_field(&m_x[0]);
	d.sort_field(&m_v[0]);
	d.sort_field(&m_a[0]);
	d.sort_field(&m_masses[0]);
	d.sort_field(&m_density[0]);
	d.sort_field(&m_particleId[0]);
	d.sort_field(&m_particleState[0]);

	if (m_viscosity)
		m_viscosity->performNeighborhoodSearchSort();
	if (m_surfaceTension)
		m_surfaceTension->performNeighborhoodSearchSort();
	if (m_vorticity)
		m_vorticity->performNeighborhoodSearchSort();
	if (m_drag)
		m_drag->performNeighborhoodSearchSort();
	if (m_elasticity)
		m_elasticity->performNeighborhoodSearchSort();
}

void SPH::FluidModel::setDensity0(const Real v)
{
	m_density0 = v; 
	initMasses(); 
}

const SPH::FieldDescription & SPH::FluidModel::getField(const std::string &name)
{
	unsigned int index = 0;
	for (auto i = 0; i < m_fields.size(); i++)
	{
		if (m_fields[i].name == name)
		{
			index = i;
			break;
		}
	}
	return m_fields[index];
}

void FluidModel::setNumActiveParticles(const unsigned int num)
{
	m_numActiveParticles = num;
}

unsigned int FluidModel::numActiveParticles() const
{
	return m_numActiveParticles;
}

void FluidModel::setDragMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_dragMethodChanged = callBackFct;
}

void FluidModel::setSurfaceMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_surfaceTensionMethodChanged = callBackFct;
}

void FluidModel::setViscosityMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_viscosityMethodChanged = callBackFct;
}

void FluidModel::setVorticityMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_vorticityMethodChanged = callBackFct;
}

void FluidModel::setElasticityMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_elasticityMethodChanged = callBackFct;
}

void FluidModel::computeSurfaceTension()
{
	if (m_surfaceTension)
		m_surfaceTension->step();
}

void FluidModel::computeViscosity()
{
	if (m_viscosity)
		m_viscosity->step();
}

void FluidModel::computeVorticity()
{
	if (m_vorticity)
		m_vorticity->step();
}

void FluidModel::computeDragForce()
{
	if (m_drag)
		m_drag->step();
}

void FluidModel::computeElasticity()
{
	if (m_elasticity)
		m_elasticity->step();
}

void FluidModel::emittedParticles(const unsigned int startIndex)
{
	if (m_viscosity)
		m_viscosity->emittedParticles(startIndex);
	if (m_surfaceTension)
		m_surfaceTension->emittedParticles(startIndex);
	if (m_vorticity)
		m_vorticity->emittedParticles(startIndex);
	if (m_drag)
		m_drag->emittedParticles(startIndex);
	if (m_elasticity)
		m_elasticity->emittedParticles(startIndex);
}

void FluidModel::setSurfaceTensionMethod(const int val)
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
		m_surfaceTension = new SurfaceTension_Becker2007(this);
	else if (m_surfaceTensionMethod == SurfaceTensionMethods::Akinci2013)
		m_surfaceTension = new SurfaceTension_Akinci2013(this);
	else if (m_surfaceTensionMethod == SurfaceTensionMethods::He2014)
		m_surfaceTension = new SurfaceTension_He2014(this);

	if (m_surfaceTension != nullptr)
		m_surfaceTension->init();

	if (m_surfaceTensionMethodChanged != nullptr)
		m_surfaceTensionMethodChanged();
}

void FluidModel::setViscosityMethod(const int val)
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
		m_viscosity = new Viscosity_Standard(this);
	else if (m_viscosityMethod == ViscosityMethods::XSPH)
		m_viscosity = new Viscosity_XSPH(this);
	else if (m_viscosityMethod == ViscosityMethods::Bender2017)
		m_viscosity = new Viscosity_Bender2017(this);
	else if (m_viscosityMethod == ViscosityMethods::Peer2015)
		m_viscosity = new Viscosity_Peer2015(this);
	else if (m_viscosityMethod == ViscosityMethods::Peer2016)
		m_viscosity = new Viscosity_Peer2016(this);
	else if (m_viscosityMethod == ViscosityMethods::Takahashi2015)
		m_viscosity = new Viscosity_Takahashi2015(this);
	else if (m_viscosityMethod == ViscosityMethods::Weiler2018)
		m_viscosity = new Viscosity_Weiler2018(this);

	if (m_viscosity != nullptr)
		m_viscosity->init();

	if (m_viscosityMethodChanged != nullptr)
		m_viscosityMethodChanged();
}


void FluidModel::setVorticityMethod(const int val)
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
		m_vorticity = new MicropolarModel_Bender2017(this);
	else if (m_vorticityMethod == VorticityMethods::VorticityConfinement)
		m_vorticity = new VorticityConfinement(this);

	if (m_vorticity != nullptr)
		m_vorticity->init();

	if (m_vorticityMethodChanged != nullptr)
		m_vorticityMethodChanged();
}

void FluidModel::setDragMethod(const int val)
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
		m_drag = new DragForce_Gissler2017(this);
	else if (m_dragMethod == DragMethods::Macklin2014)
		m_drag = new DragForce_Macklin2014(this);

	if (m_drag != nullptr)
		m_drag->init();

	if (m_dragMethodChanged != nullptr)
		m_dragMethodChanged();
}

void FluidModel::setElasticityMethod(const int val)
{
	ElasticityMethods em = static_cast<ElasticityMethods>(val);
	if ((em < ElasticityMethods::None) || (em >= ElasticityMethods::NumElasticityMethods))
		em = ElasticityMethods::None;

	if (em == m_elasticityMethod)
		return;

	delete m_elasticity;
	m_elasticity = nullptr;

	m_elasticityMethod = em;

	if (m_elasticityMethod == ElasticityMethods::Becker2009)
		m_elasticity = new Elasticity_Becker2009(this);
	else if (m_elasticityMethod == ElasticityMethods::Peer2018)
		m_elasticity = new Elasticity_Peer2018(this);

	if (m_elasticity != nullptr)
		m_elasticity->init();

	if (m_elasticityMethodChanged != nullptr)
		m_elasticityMethodChanged();
}


void FluidModel::addField(const FieldDescription &field)
{
	m_fields.push_back(field);
	std::sort(m_fields.begin(), m_fields.end(), [](FieldDescription &i, FieldDescription &j) -> bool { return (i.name < j.name); });
}

void FluidModel::removeFieldByName(const std::string &fieldName)
{
	for (auto it = m_fields.begin(); it != m_fields.end(); it++)
	{
		if (it->name == fieldName)
		{
			m_fields.erase(it);
			break;
		}
	}
}

void SPH::FluidModel::saveState(BinaryFileWriter &binWriter)
{
	binWriter.write(m_numActiveParticles);
	binWriter.writeBuffer((char*) m_particleState.data(), m_numActiveParticles * sizeof(ParticleState));

	if (m_surfaceTension)
		m_surfaceTension->saveState(binWriter);
	if (m_viscosity)
		m_viscosity->saveState(binWriter);
	if (m_vorticity)
		m_vorticity->saveState(binWriter);
	if (m_drag)
		m_drag->saveState(binWriter);
	if (m_elasticity)
		m_elasticity->saveState(binWriter);
	m_emitterSystem->saveState(binWriter);
}

void SPH::FluidModel::loadState(BinaryFileReader &binReader)
{
	binReader.read(m_numActiveParticles);
	NeighborhoodSearch *neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
	neighborhoodSearch->update_point_sets();
	neighborhoodSearch->resize_point_set(m_pointSetIndex, &getPosition(0)[0], m_numActiveParticles);

	binReader.readBuffer((char*)m_particleState.data(), m_numActiveParticles * sizeof(ParticleState));

	if (m_surfaceTension)
		m_surfaceTension->loadState(binReader);
	if (m_viscosity)
		m_viscosity->loadState(binReader);
	if (m_vorticity)
		m_vorticity->loadState(binReader);
	if (m_drag)
		m_drag->loadState(binReader);
	if (m_elasticity)
		m_elasticity->loadState(binReader);

	m_emitterSystem->loadState(binReader);
}

