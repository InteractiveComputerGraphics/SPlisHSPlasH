#include "TimeStepICSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataICSPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace std;
using namespace GenParam;

int TimeStepICSPH::LAMBDA = -1;
int TimeStepICSPH::PRESSURE_CLAMPING = -1;

TimeStepICSPH::TimeStepICSPH() :
	TimeStep()
{
	m_simulationData.init();
	m_lambda = 200000;
	m_clamping = true;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "a_ii", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getAii(fluidModelIndex, i); } });
		model->addField({ "pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); }, true });
		model->addField({ "advected density", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
		model->addField({ "pressure gradient", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureGradient(fluidModelIndex, i)[0]; } });
	}
}

TimeStepICSPH::~TimeStepICSPH(void)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("a_ii");
		model->removeFieldByName("pressure");
		model->removeFieldByName("advected density");
		model->removeFieldByName("pressure acceleration");
		model->removeFieldByName("pressure gradient");
	}
}

void TimeStepICSPH::initParameters()
{
	TimeStep::initParameters();

	LAMBDA = createNumericParameter("lambda", "Lambda", &m_lambda);
	setGroup(LAMBDA, "Simulation|ICSPH");
	setDescription(LAMBDA, "Lame parameter lambda.");
	static_cast<RealParameter*>(getParameter(LAMBDA))->setMinValue(static_cast<Real>(0.0));

	PRESSURE_CLAMPING = createBoolParameter("pressureClamping", "Enable pressure clamping", &m_clamping);
	setGroup(PRESSURE_CLAMPING, "Simulation|ICSPH");
	setDescription(PRESSURE_CLAMPING, "Turn pressure clamping on/off.");
}

void TimeStepICSPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();
	const unsigned int nModels = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// search the neighbors for all particles
	//////////////////////////////////////////////////////////////////////////
	sim->performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	//////////////////////////////////////////////////////////////////////////
	// precompute the values V_j * grad W_ij for all neighbors
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("precomputeValues")
	precomputeValues();
	STOP_TIMING_AVG
#endif

	//////////////////////////////////////////////////////////////////////////
	// compute volume/density maps boundary contribution
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	//////////////////////////////////////////////////////////////////////////
	// compute densities
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDensities(fluidModelIndex);

	//////////////////////////////////////////////////////////////////////////
	// Reset accelerations and add gravity
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		clearAccelerations(fluidModelIndex);

	//////////////////////////////////////////////////////////////////////////
	// Compute all nonpressure forces like viscosity, vorticity, ...
	//////////////////////////////////////////////////////////////////////////
	sim->computeNonPressureForces();

	//////////////////////////////////////////////////////////////////////////
	// Update the time step size, e.g. by using a CFL condition
	//////////////////////////////////////////////////////////////////////////
	sim->updateTimeStepSize();

	//////////////////////////////////////////////////////////////////////////
	// compute new velocities only considering non-pressure forces
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int m = 0; m < nModels; m++)
	{
		FluidModel *fm = sim->getFluidModel(m);
		const unsigned int numParticles = fm->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				if (fm->getParticleState(i) == ParticleState::Active)
				{
					Vector3r &vel = fm->getVelocity(i);
					vel += h * fm->getAcceleration(i);
				}
			}
		}
	}

	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		integration(fluidModelIndex);

	sim->emitParticles();
	sim->animateParticles();

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}


void TimeStepICSPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
}


/** Compute rho_adv, i ^ (0) (see equation(6) in[GHB + 20])
* using the velocities after the non-pressure forces were applied.
**/
void TimeStepICSPH::computeDensityAdv(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{

			const Real& density = model->getDensity(i);
			Real& densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
			const Vector3r& xi = model->getPosition(i);
			const Vector3r& vi = model->getVelocity(i);
			Real delta = 0.0;
			const unsigned int nFluids = sim->numberOfFluidModels();
			const unsigned int nBoundaries = sim->numberOfBoundaryModels();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Real Vj = fm_neighbor->getMass(neighborIndex) / fm_neighbor->getDensity(neighborIndex);
				const Vector3r & vj = fm_neighbor->getVelocity(neighborIndex);
				delta += Vj * (vj - vi).dot(sim->gradW(xi - xj));
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r & vj = bm_neighbor->getVelocity(neighborIndex);
					delta += bm_neighbor->getVolume(neighborIndex) * (vj - vi).dot(sim->gradW(xi - xj));
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					Vector3r vj;
					bm_neighbor->getPointVelocity(xi, vj);
					delta -= (vj - vi).dot(gradRho);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					Vector3r vj;
					bm_neighbor->getPointVelocity(xj, vj);
					delta += Vj * (vj - vi).dot(sim->gradW(xi - xj));
				);
			}

			densityAdv = density - dt * density * delta;
		}
	}
}



void TimeStepICSPH::compute_aii(const unsigned int fluidModelIndex)
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	Simulation *sim = Simulation::getCurrent();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real dt2 = dt * dt;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int) model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute a_ii
		//////////////////////////////////////////////////////////////////////////

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r &xi = model->getPosition(i);
			const Real Vi = model->getMass(i) / model->getDensity(i);
			Real Vi_Vj_gradW2 = 0.0;
			Vector3r Vj_gradW;
			Vj_gradW.setZero();
			Vector3r Vb_gradW;
			Vb_gradW.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Real Vj = fm_neighbor->getMass(neighborIndex) / fm_neighbor->getDensity(neighborIndex);
				Vj_gradW += Vj * sim->gradW(xi - xj);
				Vi_Vj_gradW2 += Vi * Vj * sim->gradW(xi - xj).squaredNorm();
			);
			
			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					Vb_gradW += bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					Vb_gradW -= gradRho;
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					Vb_gradW += Vj * sim->gradW(xi - xj);
				);
			}		

			Real& aii = m_simulationData.getAii(fluidModelIndex, i);
			aii = -density0 / m_lambda - dt2 * Vi_Vj_gradW2 - dt2 * (Vj_gradW + m_psi * Vb_gradW).dot(Vj_gradW + Vb_gradW);
		}
	}
}

void TimeStepICSPH::pressureSolve()
{
	Simulation *sim = Simulation::getCurrent();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute rho_adv
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute rho_adv,i^(0) (see equation (6) in [GHB+20])
		// using the velocities after the non-pressure forces were applied.
		//////////////////////////////////////////////////////////////////////////
		computeDensityAdv(fluidModelIndex);

		//////////////////////////////////////////////////////////////////////////
		// Compute aii (see equation (8) in [GHB+20])
		//////////////////////////////////////////////////////////////////////////
		compute_aii(fluidModelIndex);
	}


	Real avg_density_err = 0;
	m_iterations = 0;
	bool chk = false;

	while ((!chk || (m_iterations < m_minIterations)) && (m_iterations < m_maxIterations))
	{
		chk = true;
		for (unsigned int i = 0; i < nFluids; i++)
		{
			FluidModel *model = sim->getFluidModel(i);
			const Real density0 = model->getDensity0();

			avg_density_err = 0.0;
			pressureSolveIteration(i, avg_density_err);

			// Maximal allowed density fluctuation
			const Real eta = m_maxError * static_cast<Real>(0.01) * density0;  // maxError is given in percent
			chk = chk && (avg_density_err <= eta);
		}

		m_iterations++;
	}
	INCREASE_COUNTER("ICSPH - iterations", static_cast<Real>(m_iterations));
}

void TimeStepICSPH::pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	const Real density0 = model->getDensity0();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real dt2 = dt*dt;
	const Real omega = 0.5;
	Real density_error = 0.0;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			//////////////////////////////////////////////////////////////////////////
			// Compute pressure gradient (see equation (7) in [GHB+20])
			//////////////////////////////////////////////////////////////////////////
			Vector3r &grad_pi = m_simulationData.getPressureGradient(fluidModelIndex, i);
			const Real pi = m_simulationData.getPressure(fluidModelIndex, i);
			grad_pi.setZero();
			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Real Vj = fm_neighbor->getMass(neighborIndex) / fm_neighbor->getDensity(neighborIndex);
				const Real pj = m_simulationData.getPressure(pid, neighborIndex);
				grad_pi += Vj * (pi+pj) * sim->gradW(xi - xj);
			);	

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					grad_pi += m_psi * bm_neighbor->getVolume(neighborIndex) * pi * sim->gradW(xi - xj);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					grad_pi -= pi * gradRho;
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					grad_pi += Vj * pi * sim->gradW(xi - xj);
				);
			}
		}
	
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real &pi = m_simulationData.getPressure(fluidModelIndex, i);
			const Vector3r& xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Compute A p (see LHS of equation (5) in [GHB+20])
			//////////////////////////////////////////////////////////////////////////

			// compute Laplacian of p
			Real laplacian_p = 0.0;
			const Vector3r& grad_pi = m_simulationData.getPressureGradient(fluidModelIndex, i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Real Vj = fm_neighbor->getMass(neighborIndex) / fm_neighbor->getDensity(neighborIndex);
				const Vector3r &grad_pj = m_simulationData.getPressureGradient(pid, neighborIndex);
				laplacian_p += Vj * (grad_pj - grad_pi).dot(sim->gradW(xi - xj));
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					laplacian_p += bm_neighbor->getVolume(neighborIndex) * (- grad_pi).dot(sim->gradW(xi - xj));
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					laplacian_p += grad_pi.dot(gradRho);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					laplacian_p += Vj * (-grad_pi).dot(sim->gradW(xi - xj));
				);
			}

			const Real A_p = -density0 / m_lambda * pi + dt2 * laplacian_p;

			const Real &aii = m_simulationData.getAii(fluidModelIndex, i);
			const Real residuum = density0 - m_simulationData.getDensityAdv(fluidModelIndex, i) - A_p;

			if (fabs(aii) > 1.0e-5)
			{
				if (m_clamping)
					pi = std::max(pi + omega / aii * residuum, static_cast<Real>(0.0));
				else
					pi = pi + omega / aii * residuum;
			}

			//////////////////////////////////////////////////////////////////////////
			// Compute the sum of the density errors
			//////////////////////////////////////////////////////////////////////////
			density_error -= residuum;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Compute the average density error
	//////////////////////////////////////////////////////////////////////////
	avg_density_err = density_error / numParticles;
}

void TimeStepICSPH::integration(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	// Compute pressure forces
	computePressureAccels(fluidModelIndex);

	Real h = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active)
			{
				Vector3r &pos = model->getPosition(i);
				Vector3r &vel = model->getVelocity(i);
				vel += m_simulationData.getPressureAccel(fluidModelIndex, i) * h;
				pos += vel * h;
			}
		}
	}
}

void TimeStepICSPH::computePressureAccels(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = model->getDensity0();

	// Compute pressure forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = model->getPosition(i);
			const Real density_i = model->getDensity(i) / density0;

			Vector3r &ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
			ai.setZero();

			const Real density2 = density_i*density_i;
			const Real dpi = (m_simulationData.getPressure(fluidModelIndex, i)/density0) / density2;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				// Pressure 
				const Real Vj = fm_neighbor->getMass(neighborIndex) / fm_neighbor->getDensity(neighborIndex);
				const Real density_j = fm_neighbor->getDensity(neighborIndex) / fm_neighbor->getDensity0();
				const Real densityj2 = density_j*density_j;
				const Real dpj = (m_simulationData.getPressure(pid, neighborIndex)/ fm_neighbor->getDensity0())  / densityj2;
				ai -= Vj * (dpi + fm_neighbor->getDensity0() / density0 * dpj) * sim->gradW(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r a = bm_neighbor->getVolume(neighborIndex) * dpi * sim->gradW(xi - xj);
					ai -= a;
					bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					const Vector3r a = -dpi * gradRho;
					ai -= a;
					bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r a = Vj * dpi* sim->gradW(xi - xj);
					ai -= a;
					bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			}
		}
	}
}

void TimeStepICSPH::performNeighborhoodSearchSort()
{
	m_simulationData.performNeighborhoodSearchSort();
}

void TimeStepICSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepICSPH::resize()
{
	m_simulationData.init();
}
