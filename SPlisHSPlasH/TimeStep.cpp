#include "TimeStep.h"
#include "TimeManager.h"
#include "SPHKernels.h"
#include "Utilities/Timing.h"
#include "EmitterSystem.h"
#include "Simulation.h"


using namespace SPH;
using namespace std;
using namespace GenParam;

int TimeStep::SOLVER_ITERATIONS = -1;
int TimeStep::MAX_ITERATIONS = -1;
int TimeStep::MAX_ERROR = -1;


TimeStep::TimeStep()
{
	m_iterations = 0;
	m_maxIterations = 100;
	m_maxError = 0.01;
}

TimeStep::~TimeStep(void)
{
}

void TimeStep::init()
{
	initParameters();
}

void TimeStep::initParameters()
{
	ParameterObject::initParameters();

	SOLVER_ITERATIONS = createNumericParameter("iterations", "Iterations", &m_iterations);
	setGroup(SOLVER_ITERATIONS, "Simulation");
	setDescription(SOLVER_ITERATIONS, "Iterations required by the pressure solver.");
	getParameter(SOLVER_ITERATIONS)->setReadOnly(true);

	MAX_ITERATIONS = createNumericParameter("maxIterations", "Max. iterations", &m_maxIterations);
	setGroup(MAX_ITERATIONS, "Simulation");
	setDescription(MAX_ITERATIONS, "Maximal number of iterations of the pressure solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS))->setMinValue(1);

	MAX_ERROR = createNumericParameter("maxError", "Max. density error(%)", &m_maxError);
	setGroup(MAX_ERROR, "Simulation");
	setDescription(MAX_ERROR, "Maximal density error (%).");
	static_cast<RealParameter*>(getParameter(MAX_ERROR))->setMinValue(1e-6);
}

void TimeStep::clearAccelerations(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int count = model->numActiveParticles();
	const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));
	for (unsigned int i=0; i < count; i++)
	{
		// Clear accelerations of dynamic particles
		if (model->getMass(i) != 0.0)
		{
			Vector3r &a = model->getAcceleration(i);
			a = grav;
		}
	}
}

void TimeStep::computeDensities(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			Real &density = model->getDensity(i);

			// Compute current density for particle i
			density = model->getVolume(i) * sim->W_zero();
			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				density += fm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			forall_boundary_neighbors(				
				// Boundary: Akinci2012
				density += bm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			)

			density *= density0;
		}
	}
}

void TimeStep::reset()
{
	m_iterations = 0;
}
