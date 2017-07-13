#include "TimeStep.h"
#include "TimeManager.h"
#include "SPHKernels.h"
#include "PBF/SimulationDataPBF.h"
#include "SPlisHSPlasH/Utilities/Timing.h"
#include "SurfaceTension/SurfaceTension_Becker2007.h"
#include "SurfaceTension/SurfaceTension_Akinci2013.h"
#include "SurfaceTension/SurfaceTension_He2014.h"
#include "Viscosity/Viscosity_XSPH.h"
#include "Viscosity/Viscosity_Standard.h"
#include "Viscosity/Viscosity_Bender2017.h"
#include "Vorticity/VorticityConfinement.h"
#include "EmitterSystem.h"


using namespace SPH;
using namespace std;

TimeStep::TimeStep(FluidModel *model)
{
	m_model = model;
	m_iterations = 0;
	m_iterationsV = 0;
	m_cflMethod = 1;
	m_cflFactor = 0.5;
	m_cflMaxTimeStepSize = 0.005;
	m_maxIterations = 100;
	m_maxError = 0.01;
	m_maxIterationsV = 100;
	m_maxErrorV = 0.1;
	m_viscosity = NULL;
	m_viscosityMethod = ViscosityMethods::None;
	setViscosityMethod(ViscosityMethods::XSPH);
	m_surfaceTension = NULL;
	m_surfaceTensionMethod = SurfaceTensionMethods::None;
	m_fluidModel = VorticityMethods::None;
	m_vorticity = NULL;
}

TimeStep::~TimeStep(void)
{
}

void TimeStep::clearAccelerations()
{
	const unsigned int count = m_model->numActiveParticles();
	const Vector3r &grav = m_model->getGravitation();
	for (unsigned int i=0; i < count; i++)
	{
		// Clear accelerations of dynamic particles
		if (m_model->getMass(i) != 0.0)
		{
			Vector3r &a = m_model->getAcceleration(i);
			a = grav;
		}
	}
}

void TimeStep::computeDensities()
{
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real density0 = m_model->getDensity0();
	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			Real &density = m_model->getDensity(i);

			// Compute current density for particle i
			density = m_model->getMass(i) * m_model->W_zero();
			const Vector3r &xi = m_model->getPosition(0, i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				density += m_model->getMass(neighborIndex) * m_model->W(xi - xj);
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int pid=1; pid < m_model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
					const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
					
					// Boundary: Akinci2012
					density += m_model->getBoundaryPsi(pid, neighborIndex) * m_model->W(xi - xj);
				}
			}
		}
	}
}

void TimeStep::updateTimeStepSize()
{
	if (m_cflMethod == 1)
		updateTimeStepSizeCFL(0.0001);
	else if(m_cflMethod == 2)
	{
		Real h = TimeManager::getCurrent()->getTimeStepSize();
		updateTimeStepSizeCFL(0.0001);
		if (m_iterations > 10)
			h *= 0.9;
		else if (m_iterations < 5)
			h *= 1.1;
		h = min(h, TimeManager::getCurrent()->getTimeStepSize());
		TimeManager::getCurrent()->setTimeStepSize(h);
	}
}

void TimeStep::updateTimeStepSizeCFL(const Real minTimeStepSize)
{
	const Real radius = m_model->getParticleRadius();
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

void TimeStep::computeNonPressureForces()
{
	computeSurfaceTension();
	computeViscosity();
	computeVorticity();
}

void TimeStep::computeSurfaceTension()
{
	if (m_surfaceTension)
		m_surfaceTension->step();
}

void TimeStep::computeViscosity()
{
	if (m_viscosity)
		m_viscosity->step();
}

void TimeStep::computeVorticity()
{
	if (m_vorticity)
		m_vorticity->step();
}

void TimeStep::performNeighborhoodSearch()
{
	START_TIMING("neighborhood_search");
	m_model->getNeighborhoodSearch()->find_neighbors();
	STOP_TIMING_AVG;
}

void TimeStep::reset()
{
	m_model->reset();
	if (m_surfaceTension)
		m_surfaceTension->reset();
	if (m_viscosity)
		m_viscosity->reset();
	m_iterations = 0;

	TimeManager::getCurrent()->setTimeStepSize(0.001);
}

void TimeStep::setSurfaceTensionMethod(SurfaceTensionMethods val)
{
	if ((val < SurfaceTensionMethods::None) || (val > SurfaceTensionMethods::He2014))
		val = SurfaceTensionMethods::None;

	if (val == m_surfaceTensionMethod)
		return;

	delete m_surfaceTension;
	m_surfaceTension = NULL;

	m_surfaceTensionMethod = val;
	if (m_surfaceTensionMethod == SurfaceTensionMethods::Becker2007)
		m_surfaceTension = new SurfaceTension_Becker2007(m_model);
	else if (m_surfaceTensionMethod == SurfaceTensionMethods::Akinci2013)
		m_surfaceTension = new SurfaceTension_Akinci2013(m_model);
	else if (m_surfaceTensionMethod == SurfaceTensionMethods::He2014)
		m_surfaceTension = new SurfaceTension_He2014(m_model);
}

void SPH::TimeStep::setViscosityMethod(ViscosityMethods val)
{
	if ((val < ViscosityMethods::None) || (val > ViscosityMethods::Bender2017))
		val = ViscosityMethods::XSPH;

	if (val == m_viscosityMethod)
		return;

	delete m_viscosity;
	m_viscosity = NULL;

	m_viscosityMethod = val;

	if (m_viscosityMethod == ViscosityMethods::Standard)
		m_viscosity = new Viscosity_Standard(m_model);	
	else if (m_viscosityMethod == ViscosityMethods::XSPH)
		m_viscosity = new Viscosity_XSPH(m_model);
	else if (m_viscosityMethod == ViscosityMethods::Bender2017)
		m_viscosity = new Viscosity_Bender2017(m_model);
}

void TimeStep::setVorticityMethod(SPH::VorticityMethods val)
{
	if ((val < VorticityMethods::None) || (val > VorticityMethods::VorticityConfinement))
		val = VorticityMethods::None;

	if (val == m_fluidModel)
		return;

	delete m_vorticity;
	m_vorticity = NULL;

	m_fluidModel = val;

	if (m_fluidModel == VorticityMethods::Micropolar)
		m_vorticity = new MicropolarModel_Bender2017(m_model);
	else if (m_fluidModel == VorticityMethods::VorticityConfinement)
		m_vorticity = new VorticityConfinement(m_model);
}

void TimeStep::emitParticles()
{
	getModel()->getEmitterSystem().step(this);
}
