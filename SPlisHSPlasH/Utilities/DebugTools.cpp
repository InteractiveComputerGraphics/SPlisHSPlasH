#include "DebugTools.h"

#include "Utilities/Logger.h"
#include <SPlisHSPlasH/Simulation.h>

using namespace SPH;
using namespace GenParam;

int DebugTools::DETERMINE_THREAD_IDS = -1;
int DebugTools::DETERMINE_NUM_NEIGHBORS = -1;
int DebugTools::DETERMINE_VELOCITY_CHANGES = -1;

DebugTools::DebugTools() :
	ParameterObject()
{
	m_determineThreadIds = false;

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "threadId", FieldType::UInt, [this, fluidModelIndex](const unsigned int i) -> unsigned int* { return &m_threadIds[fluidModelIndex][i]; } });
		model->addField({ "numNeighbors", FieldType::UInt, [this, fluidModelIndex](const unsigned int i) -> unsigned int* { return &m_numNeighbors[fluidModelIndex][i]; } });
		model->addField({ "velocityChanges", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_velocityChanges[fluidModelIndex][i][0]; } });
	}
}

DebugTools::~DebugTools()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("threadId");
		model->removeFieldByName("numNeighbors");
		model->removeFieldByName("velocityChanges");
	}
}

void DebugTools::init()
{
	initParameters();

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_threadIds.resize(nModels);
	m_numNeighbors.resize(nModels);
	m_vOld.resize(nModels);
	m_velocityChanges.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		m_threadIds[i].resize(fm->numParticles(), 0);
		m_numNeighbors[i].resize(fm->numParticles(), 0);
		m_vOld[i].resize(fm->numParticles(), Vector3r::Zero());
		m_velocityChanges[i].resize(fm->numParticles(), Vector3r::Zero());
	}
}

void SPH::DebugTools::cleanup()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_threadIds[i].clear();
		m_numNeighbors[i].clear();
		m_vOld[i].clear();
		m_velocityChanges[i].clear();
	}
	m_threadIds.clear();
	m_numNeighbors.clear();
	m_vOld.clear();
	m_velocityChanges.clear();
}

void DebugTools::initParameters()
{
	ParameterObject::initParameters();

	DETERMINE_THREAD_IDS = createBoolParameter("determineThreadIds", "Determine Thread IDs", &m_determineThreadIds);
	setGroup(DETERMINE_THREAD_IDS, "Debug Tools");
	setDescription(DETERMINE_THREAD_IDS, "Determine Thread IDs and add a corresponding particle field.");

	DETERMINE_NUM_NEIGHBORS = createBoolParameter("determineNumNeighbors", "Determine # neighbors", &m_determineNumNeighbors);
	setGroup(DETERMINE_NUM_NEIGHBORS, "Debug Tools");
	setDescription(DETERMINE_NUM_NEIGHBORS, "Determine number of neighbors and add a corresponding particle field.");

	DETERMINE_VELOCITY_CHANGES = createBoolParameter("determineVelocityChanges", "Determine velocity changes", &m_determineVelocityChanges);
	setGroup(DETERMINE_VELOCITY_CHANGES, "Debug Tools");
	setDescription(DETERMINE_VELOCITY_CHANGES, "Determine velocity change of each particle and add a corresponding particle field.");
}

void SPH::DebugTools::determineThreadIds()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
#ifdef _OPENMP
				m_threadIds[fluidModelIndex][i] = omp_get_thread_num();
#else
				m_threadIds[fluidModelIndex][i] = 0;
#endif
			}
		}
	}
}

void SPH::DebugTools::determineNumNeighbors()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				m_numNeighbors[fluidModelIndex][i] = 0;
				for (unsigned int pid = 0; pid < sim->numberOfPointSets(); pid++)
					m_numNeighbors[fluidModelIndex][i] += sim->numberOfNeighbors(fluidModelIndex, pid, i);
			}
		}
	}
}

void SPH::DebugTools::determineVelocityChanges()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				m_velocityChanges[fluidModelIndex][i] = model->getVelocity(i) - m_vOld[fluidModelIndex][i];
				m_vOld[fluidModelIndex][i] = model->getVelocity(i);
			}
		}
	}
}

void DebugTools::step()
{
	if (m_determineThreadIds)
		determineThreadIds();

	if (m_determineNumNeighbors)
		determineNumNeighbors();

	if (m_determineVelocityChanges)
		determineVelocityChanges();
}

void DebugTools::reset()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numParticles(); j++)
		{
			m_threadIds[i][j] = 0;
			m_numNeighbors[i][j] = 0;
			m_vOld[i][j].setZero();
			m_velocityChanges[i][j].setZero();
		}
	}
}

void DebugTools::performNeighborhoodSearchSort()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		const unsigned int numPart = fm->numActiveParticles();
		if (numPart != 0)
		{
			auto const& d = sim->getNeighborhoodSearch()->point_set(fm->getPointSetIndex());
			d.sort_field(&m_vOld[i][0]);
		}
	}
}

void DebugTools::emittedParticles(FluidModel* model, const unsigned int startIndex)
{
}

