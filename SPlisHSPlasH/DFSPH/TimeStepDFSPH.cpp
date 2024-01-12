#include "TimeStepDFSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataDFSPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace std;
using namespace GenParam;


int TimeStepDFSPH::SOLVER_ITERATIONS_V = -1;
int TimeStepDFSPH::MAX_ITERATIONS_V = -1;
int TimeStepDFSPH::MAX_ERROR_V = -1;
int TimeStepDFSPH::USE_DIVERGENCE_SOLVER = -1;


TimeStepDFSPH::TimeStepDFSPH() :
	TimeStep(),
	m_simulationData()
{
	m_simulationData.init();
	m_iterationsV = 0;
	m_enableDivergenceSolver = true;
	m_maxIterationsV = 100;
	m_maxErrorV = static_cast<Real>(0.1);

	// add particle fields - then they can be used for the visualization and export
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "factor", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getFactor(fluidModelIndex, i); } });
		model->addField({ "advected density", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "p / rho^2", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureRho2(fluidModelIndex, i); }, true });
		model->addField({ "p_v / rho^2", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureRho2_V(fluidModelIndex, i); }, true });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}
}

TimeStepDFSPH::~TimeStepDFSPH(void)
{
	// remove all particle fields
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("factor");
		model->removeFieldByName("advected density");
		model->removeFieldByName("p / rho^2");
		model->removeFieldByName("p_v / rho^2");
		model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepDFSPH::initParameters()
{
	TimeStep::initParameters();

	SOLVER_ITERATIONS_V = createNumericParameter("iterationsV", "Iterations (divergence)", &m_iterationsV);
	setGroup(SOLVER_ITERATIONS_V, "Simulation|DFSPH");
	setDescription(SOLVER_ITERATIONS_V, "Iterations required by the divergence solver.");
	getParameter(SOLVER_ITERATIONS_V)->setReadOnly(true);

	MAX_ITERATIONS_V = createNumericParameter("maxIterationsV", "Max. iterations (divergence)", &m_maxIterationsV);
	setGroup(MAX_ITERATIONS_V, "Simulation|DFSPH");
	setDescription(MAX_ITERATIONS_V, "Maximal number of iterations of the divergence solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS_V))->setMinValue(1);

	MAX_ERROR_V = createNumericParameter("maxErrorV", "Max. divergence error(%)", &m_maxErrorV);
	setGroup(MAX_ERROR_V, "Simulation|DFSPH");
	setDescription(MAX_ERROR_V, "Maximal divergence error (%).");
	static_cast<RealParameter*>(getParameter(MAX_ERROR_V))->setMinValue(static_cast<Real>(1e-6));

	USE_DIVERGENCE_SOLVER = createBoolParameter("enableDivergenceSolver", "Enable divergence solver", &m_enableDivergenceSolver);
	setGroup(USE_DIVERGENCE_SOLVER, "Simulation|DFSPH");
	setDescription(USE_DIVERGENCE_SOLVER, "Turn divergence solver on/off.");
}

void TimeStepDFSPH::step()
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
	// Compute the factor alpha_i for all particles i
	// using the equation (11) in [BK17]
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("computeDFSPHFactor");
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDFSPHFactor(fluidModelIndex);
	STOP_TIMING_AVG;

	//////////////////////////////////////////////////////////////////////////
	// Perform divergence solve (see Algorithm 2 in [BK17])
	//////////////////////////////////////////////////////////////////////////
	if (m_enableDivergenceSolver)
	{
		START_TIMING("divergenceSolve");
		divergenceSolve();
		STOP_TIMING_AVG
	}
	else
		m_iterationsV = 0;

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

	//////////////////////////////////////////////////////////////////////////
	// Perform constant density solve (see Algorithm 3 in [BK17])
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	//////////////////////////////////////////////////////////////////////////
	// compute final positions
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
					Vector3r &xi = fm->getPosition(i);
					const Vector3r &vi = fm->getVelocity(i);
					xi += h * vi;
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// emit new particles and perform an animation field step
	//////////////////////////////////////////////////////////////////////////
	sim->emitParticles();
	sim->animateParticles();

	//////////////////////////////////////////////////////////////////////////
	// Compute new time
	//////////////////////////////////////////////////////////////////////////
	tm->setTime (tm->getTime () + h);
}


void TimeStepDFSPH::pressureSolve()
{
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real invH = static_cast<Real>(1.0) / h;
	const Real invH2 = static_cast<Real>(1.0) / h2;
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute rho_adv
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const Real density0 = model->getDensity0();
		const int numParticles = (int)model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Compute rho_adv,i^(0) (see equation in Section 3.3 in [BK17])
				// using the velocities after the non-pressure forces were applied.
				//////////////////////////////////////////////////////////////////////////
				computeDensityAdv(fluidModelIndex, i, h, density0);

				//////////////////////////////////////////////////////////////////////////
				// In the end of Section 3.3 [BK17] we have to multiply the density 
				// error with the factor alpha_i divided by h^2. Hence, we multiply 
				// the factor directly by 1/h^2 here.
				//////////////////////////////////////////////////////////////////////////
				m_simulationData.getFactor(fluidModelIndex, i) *= invH2;

				//////////////////////////////////////////////////////////////////////////
				// For the warm start we use 0.5 times the old pressure value.
				// Note: We divide the value by h^2 since we multiplied it by h^2 at the end of 
				// the last time step to make it independent of the time step size.
				//////////////////////////////////////////////////////////////////////////
#ifdef USE_WARMSTART
				if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 1.0)
					m_simulationData.getPressureRho2(fluidModelIndex, i) = static_cast<Real>(0.5) * min(m_simulationData.getPressureRho2(fluidModelIndex, i), static_cast<Real>(0.00025)) * invH2;
				else 
					m_simulationData.getPressureRho2(fluidModelIndex, i) = 0.0;
#else 
				//////////////////////////////////////////////////////////////////////////
				// If we don't use a warm start, we directly compute a pressure value
				// by multiplying the density error with the factor.
				//////////////////////////////////////////////////////////////////////////
				//m_simulationData.getPressureRho2(fluidModelIndex, i) = 0.0;
				const Real s_i = static_cast<Real>(1.0) - m_simulationData.getDensityAdv(fluidModelIndex, i);
				const Real residuum = min(s_i, static_cast<Real>(0.0));     // r = b - A*p
				m_simulationData.getPressureRho2(fluidModelIndex, i) = -residuum * m_simulationData.getFactor(fluidModelIndex, i);
#endif
			}
		}
	}

	m_iterations = 0;

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////
	
	Real avg_density_err = 0.0;
	bool chk = false;


	//////////////////////////////////////////////////////////////////////////
	// Perform solver iterations
	//////////////////////////////////////////////////////////////////////////
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

	INCREASE_COUNTER("DFSPH - iterations", static_cast<Real>(m_iterations));

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();
		const Real density0 = model->getDensity0();
		
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Time integration of the pressure accelerations to get new velocities
				//////////////////////////////////////////////////////////////////////////
				computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2Data(), true);
				model->getVelocity(i) += h * m_simulationData.getPressureAccel(fluidModelIndex, i);
			}
		}
	}
#ifdef USE_WARMSTART
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Multiply by h^2, the time step size has to be removed 
				// to make the pressure value independent 
				// of the time step size
				//////////////////////////////////////////////////////////////////////////
				m_simulationData.getPressureRho2(fluidModelIndex, i) *= h2;
			}		
		}
	}
#endif
}

void TimeStepDFSPH::divergenceSolve()
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	Simulation *sim = Simulation::getCurrent();
	const unsigned int maxIter = m_maxIterationsV;
	const Real maxError = m_maxErrorV;
	const unsigned int nFluids = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute divergence of velocity field
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();

		#pragma omp parallel default(shared)
		{		
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Compute rho_adv,i^(0) (see equation (9) in Section 3.2 [BK17])
				// using the velocities after the non-pressure forces were applied.
				//////////////////////////////////////////////////////////////////////////
				computeDensityChange(fluidModelIndex, i, h);

				Real densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
				densityAdv = max(densityAdv, static_cast<Real>(0.0));

				unsigned int numNeighbors = 0;
				for (unsigned int pid = 0; pid < sim->numberOfPointSets(); pid++)
					numNeighbors += sim->numberOfNeighbors(fluidModelIndex, pid, i);

				// in case of particle deficiency do not perform a divergence solve
				if (!sim->is2DSimulation())
				{
					if (numNeighbors < 20)
						densityAdv = 0.0;
				}
				else
				{
					if (numNeighbors < 7)
						densityAdv = 0.0;
				}
				
				//////////////////////////////////////////////////////////////////////////
				// In equation (11) [BK17] we have to multiply the divergence 
				// error with the factor divided by h. Hence, we multiply the factor
				// directly by 1/h here.
				//////////////////////////////////////////////////////////////////////////
				m_simulationData.getFactor(fluidModelIndex, i) *= invH;

				//////////////////////////////////////////////////////////////////////////
				// For the warm start we use 0.5 times the old pressure value.
				// Divide the value by h. We multiplied it by h at the end of 
				// the last time step to make it independent of the time step size.
				//////////////////////////////////////////////////////////////////////////
#ifdef USE_WARMSTART_V
				if (densityAdv > 0.0)
					m_simulationData.getPressureRho2_V(fluidModelIndex, i) = static_cast<Real>(0.5) * min(m_simulationData.getPressureRho2_V(fluidModelIndex, i), static_cast<Real>(0.5)) * invH;
				else
					m_simulationData.getPressureRho2_V(fluidModelIndex, i) = 0.0;
#else 
				//////////////////////////////////////////////////////////////////////////
				// If we don't use a warm start, directly compute a pressure value
				// by multiplying the divergence error with the factor.
				//////////////////////////////////////////////////////////////////////////
				m_simulationData.getPressureRho2_V(fluidModelIndex, i) = densityAdv * m_simulationData.getFactor(fluidModelIndex, i);
#endif
			}
		}
	}

	m_iterationsV = 0;

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////
	
	Real avg_density_err = 0.0;
	bool chk = false;

	//////////////////////////////////////////////////////////////////////////
	// Perform solver iterations
	//////////////////////////////////////////////////////////////////////////
	while ((!chk || (m_iterationsV < 1)) && (m_iterationsV < maxIter))
	{
		chk = true;
		for (unsigned int i = 0; i < nFluids; i++)
		{
			FluidModel *model = sim->getFluidModel(i);
			const Real density0 = model->getDensity0();

			avg_density_err = 0.0;
			divergenceSolveIteration(i, avg_density_err);

			// Maximal allowed density fluctuation
			// use maximal density error divided by time step size
			const Real eta = (static_cast<Real>(1.0) / h) * maxError * static_cast<Real>(0.01) * density0;  // maxError is given in percent
			chk = chk && (avg_density_err <= eta);
		}

		m_iterationsV++;
	}

	INCREASE_COUNTER("DFSPH - iterationsV", static_cast<Real>(m_iterationsV));

	
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();
		const Real density0 = model->getDensity0();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Time integration of the pressure accelerations
				//////////////////////////////////////////////////////////////////////////
				computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2VData(), true);
				model->getVelocity(i) += h * m_simulationData.getPressureAccel(fluidModelIndex, i);

				m_simulationData.getFactor(fluidModelIndex, i) *= h;
			}
		}
	}
#ifdef USE_WARMSTART_V
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Multiply by h, the time step size has to be removed 
				// to make the pressure value independent 
				// of the time step size
				//////////////////////////////////////////////////////////////////////////		
				m_simulationData.getPressureRho2_V(fluidModelIndex, i) *= h;
			}
		}
	}
#endif
}


void TimeStepDFSPH::pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	
	Real density_error = 0.0;

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure accelerations using the current pressure values.
		// (see Algorithm 3, line 7 in [BK17])
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2Data());
		}

		//////////////////////////////////////////////////////////////////////////
		// Update pressure values
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{						
			if (model->getParticleState(i) != ParticleState::Active)
				continue;
				
			Real aij_pj = compute_aij_pj(fluidModelIndex, i);
			aij_pj *= h * h;

			//////////////////////////////////////////////////////////////////////////
			// Compute source term: s_i = 1 - rho_adv
			// Note: that due to our multiphase handling, the multiplier rho0
			// is missing here
			//////////////////////////////////////////////////////////////////////////
			const Real& densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
			const Real s_i = static_cast<Real>(1.0) - densityAdv;


			//////////////////////////////////////////////////////////////////////////
			// Update the value p/rho^2 (in [BK17] this is kappa/rho):
			// 
			// alpha_i = -1 / (a_ii * rho_i^2)
			// p_rho2_i = (p_i / rho_i^2)
			// 
			// Therefore, the following lines compute the Jacobi iteration:
			// p_i := p_i + 1/a_ii (source_term_i - a_ij * p_j)
			//////////////////////////////////////////////////////////////////////////
			Real& p_rho2_i = m_simulationData.getPressureRho2(fluidModelIndex, i);
			const Real residuum = min(s_i - aij_pj, static_cast<Real>(0.0));     // r = b - A*p
			//p_rho2_i -= residuum * m_simulationData.getFactor(fluidModelIndex, i);

			p_rho2_i = max(p_rho2_i - 0.5 * (s_i - aij_pj) * m_simulationData.getFactor(fluidModelIndex, i), 0.0);

			//////////////////////////////////////////////////////////////////////////
			// Compute the sum of the density errors
			//////////////////////////////////////////////////////////////////////////
			density_error -= density0 * residuum;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Compute the average density error
	//////////////////////////////////////////////////////////////////////////
	avg_density_err = density_error / numParticles;
}

void TimeStepDFSPH::divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	
	Real density_error = 0.0;
	
	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure accelerations using the current pressure values.
 		// (see Algorithm 2, line 7 in [BK17])
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2VData());
		}

		//////////////////////////////////////////////////////////////////////////
		// Update pressure 
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			Real aij_pj = compute_aij_pj(fluidModelIndex, i);
			aij_pj *= h;

			//////////////////////////////////////////////////////////////////////////
			// Compute source term: s_i = -d rho / dt
			//////////////////////////////////////////////////////////////////////////
			const Real& densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
			const Real s_i = -densityAdv;

			//////////////////////////////////////////////////////////////////////////
			// Update the value p/rho^2:
			// 
			// alpha_i = -1 / (a_ii * rho_i^2)
			// pv_rho2_i = (pv_i / rho_i^2)
			// 
			// Therefore, the following line computes the Jacobi iteration:
			// pv_i := pv_i + 1/a_ii (source_term_i - a_ij * pv_j)
			//////////////////////////////////////////////////////////////////////////
			Real& pv_rho2_i = m_simulationData.getPressureRho2_V(fluidModelIndex, i);
			Real residuum = min(s_i - aij_pj, static_cast<Real>(0.0));     // r = b - A*p

			unsigned int numNeighbors = 0;
			for (unsigned int pid = 0; pid < sim->numberOfPointSets(); pid++)
				numNeighbors += sim->numberOfNeighbors(fluidModelIndex, pid, i);

			// in case of particle deficiency do not perform a divergence solve
			if (!sim->is2DSimulation())
			{
				if (numNeighbors < 20)
					residuum = 0.0;
			}
			else
			{
				if (numNeighbors < 7)
					residuum = 0.0;
			}
			//pv_rho2_i -= residuum * m_simulationData.getFactor(fluidModelIndex, i);
			pv_rho2_i = max(pv_rho2_i - 0.5*(s_i - aij_pj) * m_simulationData.getFactor(fluidModelIndex, i), 0.0);


			//////////////////////////////////////////////////////////////////////////
			// Compute the sum of the divergence errors
			//////////////////////////////////////////////////////////////////////////
			density_error -= density0 * residuum;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Compute the average divergence error
	//////////////////////////////////////////////////////////////////////////
	avg_density_err = density_error / numParticles;
}



void TimeStepDFSPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_iterations = 0;
	m_iterationsV = 0;
}

void TimeStepDFSPH::performNeighborhoodSearchSort()
{
	m_simulationData.performNeighborhoodSearchSort();
}

void TimeStepDFSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepDFSPH::resize()
{
	m_simulationData.init();
}

#ifdef USE_AVX

void TimeStepDFSPH::computeDFSPHFactor(const unsigned int fluidModelIndex)
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const int numParticles = (int)model->numActiveParticles();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure stiffness denominator
		//////////////////////////////////////////////////////////////////////////

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			//////////////////////////////////////////////////////////////////////////
			// Compute gradient dp_i/dx_j * (1/kappa)  and dp_j/dx_j * (1/kappa)
			// (see Equation (8) and the previous one [BK17])
			// Note: That in all quantities rho0 is missing due to our
			// implementation of multiphase simulations.
			//////////////////////////////////////////////////////////////////////////
			const Vector3r& xi = model->getPosition(i);

			Real sum_grad_p_k;
			Vector3r grad_p_i;
			Vector3f8 xi_avx(xi);
			Scalarf8 sum_grad_p_k_avx(0.0f);
			Vector3f8 grad_p_i_avx;
			grad_p_i_avx.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_avx_nox(
				compute_xj(fm_neighbor, pid);
				compute_Vj(fm_neighbor);
				compute_Vj_gradW();
				sum_grad_p_k_avx += V_gradW.squaredNorm();
				grad_p_i_avx += V_gradW;
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors_avx(
					const Scalarf8 V_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
					const Vector3f8 grad_p_j = CubicKernel_AVX::gradW(xj_avx - xi_avx) * V_avx;
					grad_p_i_avx -= grad_p_j;
				);
			}

			sum_grad_p_k = sum_grad_p_k_avx.reduce();
			grad_p_i[0] = grad_p_i_avx.x().reduce();
			grad_p_i[1] = grad_p_i_avx.y().reduce();
			grad_p_i[2] = grad_p_i_avx.z().reduce();

			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					grad_p_i -= gradRho;
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
					grad_p_i -= grad_p_j;
				);
			}

			sum_grad_p_k += grad_p_i.squaredNorm();

			//////////////////////////////////////////////////////////////////////////
			// Compute factor alpha_i / rho_i (see Equation (11) in [BK17])
			//////////////////////////////////////////////////////////////////////////
			Real& factor = m_simulationData.getFactor(fluidModelIndex, i);
			if (sum_grad_p_k > m_eps)
				factor = static_cast<Real>(1.0) / (sum_grad_p_k);
			else
				factor = 0.0;
		}
	}
}

/** Compute rho_adv,i^(0) (see equation in Section 3.3 in [BK17])
  * using the velocities after the non-pressure forces were applied.
**/
void TimeStepDFSPH::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const Real h, const Real density0)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real &density = model->getDensity(i);
	Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	const Vector3r &xi = model->getPosition(i);
	const Vector3r &vi = model->getVelocity(i);
	Real delta = 0.0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	Scalarf8 delta_avx(0.0f);
	const Vector3f8 xi_avx(xi);
	Vector3f8 vi_avx(vi);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_avx_nox(
		compute_xj(fm_neighbor, pid);
		compute_Vj(fm_neighbor);
		compute_Vj_gradW();
		const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getVelocity(0), count);
		delta_avx += (vi_avx - vj_avx).dot(V_gradW);
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors_avx(
			const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
			const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVelocity(0), count);
			delta_avx += Vj_avx * (vi_avx - vj_avx).dot(CubicKernel_AVX::gradW(xi_avx - xj_avx));
		);
	}

	delta = delta_avx.reduce();

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xi, vj);
			delta -= (vi - vj).dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xj, vj);
			delta += Vj * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}

	densityAdv = density / density0 + h*delta;
}

/** Compute rho_adv,i^(0) (see equation (9) in Section 3.2 [BK17])
  * using the velocities after the non-pressure forces were applied.
  */
void TimeStepDFSPH::computeDensityChange(const unsigned int fluidModelIndex, const unsigned int i, const Real h)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Vector3r &xi = model->getPosition(i);
	const Vector3r &vi = model->getVelocity(i);
	unsigned int numNeighbors = 0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	Scalarf8 densityAdv_avx(0.0f);
	const Vector3f8 xi_avx(xi);
	Vector3f8 vi_avx(vi);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_avx_nox(
		compute_xj(fm_neighbor, pid);
		compute_Vj(fm_neighbor);
		compute_Vj_gradW();
		const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getVelocity(0), count);
		densityAdv_avx += (vi_avx - vj_avx).dot(V_gradW);
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors_avx(
			const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
			const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVelocity(0), count);
			densityAdv_avx += Vj_avx * (vi_avx - vj_avx).dot(CubicKernel_AVX::gradW(xi_avx - xj_avx));
		);
	}

	Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	densityAdv = densityAdv_avx.reduce();

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xi, vj);
			densityAdv -= (vi - vj).dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xj, vj);
			densityAdv += Vj * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}
}

/** Compute pressure accelerations using the current pressure values of the particles
 */
void TimeStepDFSPH::computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, std::vector<std::vector<Real>>& pressure_rho2, const bool applyBoundaryForces)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);

	if (model->getParticleState(i) != ParticleState::Active)
		return;

	// p_rho2_i = (p_i / rho_i^2)
	const Real p_rho2_i = pressure_rho2[fluidModelIndex][i];
	const Vector3r &xi = model->getPosition(i);

	Scalarf8 p_rho2_i_avx(p_rho2_i);
	const Vector3f8 xi_avx(xi);
	Vector3f8 delta_ai_avx;
	delta_ai_avx.setZero();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_avx_nox(
		compute_xj(fm_neighbor, pid);
		compute_Vj(fm_neighbor);
		compute_Vj_gradW();
		const Scalarf8 densityFrac_avx(fm_neighbor->getDensity0() / density0);

		// p_rho2_j = (p_j / rho_j^2)
		const Scalarf8 p_rho2_j_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &pressure_rho2[pid][0], count);
		const Scalarf8 pSum = p_rho2_i_avx + densityFrac_avx * p_rho2_j_avx;
		delta_ai_avx -= V_gradW * pSum;
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (fabs(p_rho2_i) > m_eps)
	{
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			const Scalarf8 mi_avx(model->getMass(i));
			forall_boundary_neighbors_avx(
				const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);

				// Directly update velocities instead of storing pressure accelerations
				const Vector3f8 a = -CubicKernel_AVX::gradW(xi_avx - xj_avx) * (Vj_avx * p_rho2_i_avx);
				delta_ai_avx += a;

				if (applyBoundaryForces)
					bm_neighbor->addForce(xj_avx, -a * mi_avx, count);
			);
		}
	}

	ai[0] = delta_ai_avx.x().reduce();
	ai[1] = delta_ai_avx.y().reduce();
	ai[2] = delta_ai_avx.z().reduce();

	if (fabs(p_rho2_i) > m_eps)
	{
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		{
			forall_density_maps(
				const Vector3r a = (Real) 1.0 * p_rho2_i * gradRho;			
				ai += a;

				if (applyBoundaryForces)
					bm_neighbor->addForce(xj, -model->getMass(i) * a);
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		{
			forall_volume_maps(
				const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
				const Vector3r a = (Real) 1.0 * p_rho2_i * grad_p_j;			
				ai += a;

				if (applyBoundaryForces)
					bm_neighbor->addForce(xj, -model->getMass(i) * a);
			);
		}
	}
}


Real TimeStepDFSPH::compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute A*p which is the change of the density when applying the 
	// pressure forces. 
	// \sum_j a_ij * p_j = h^2 \sum_j V_j (a_i - a_j) * gradW_ij
	// This is the RHS of Equation (12) in [BK17]
	//////////////////////////////////////////////////////////////////////////
	const Vector3r& xi = model->getPosition(i);
	const Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
	const Vector3f8 xi_avx(xi);
	const Vector3f8 ai_avx(ai);
	Scalarf8 aij_pj_avx;
	aij_pj_avx.setZero();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_avx_nox(
		compute_xj(fm_neighbor, pid);
		compute_Vj(fm_neighbor);
		compute_Vj_gradW();

		const Vector3f8 aj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &m_simulationData.getPressureAccel(pid, 0), count);
		aij_pj_avx += (ai_avx - aj_avx).dot(V_gradW);
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors_avx(
			const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
			aij_pj_avx += Vj_avx * ai_avx.dot(CubicKernel_AVX::gradW(xi_avx - xj_avx));
		);
	}

	Real aij_pj = aij_pj_avx.reduce();

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			aij_pj -= ai.dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			aij_pj += Vj * ai.dot(sim->gradW(xi - xj));
		);
	}
	return aij_pj;
}



#else

void TimeStepDFSPH::computeDFSPHFactor(const unsigned int fluidModelIndex)
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const int numParticles = (int) model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure stiffness denominator
		//////////////////////////////////////////////////////////////////////////

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			//////////////////////////////////////////////////////////////////////////
			// Compute gradient dp_i/dx_j * (1/kappa)  and dp_j/dx_j * (1/kappa)
			// (see Equation (8) and the previous one [BK17])
			// Note: That in all quantities rho0 is missing due to our
			// implementation of multiphase simulations.
			//////////////////////////////////////////////////////////////////////////
			const Vector3r &xi = model->getPosition(i);
			Real sum_grad_p_k = 0.0;
			Vector3r grad_p_i;
			grad_p_i.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r grad_p_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
				sum_grad_p_k += grad_p_j.squaredNorm();
				grad_p_i -= grad_p_j;
			);
			
			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
					grad_p_i -= grad_p_j;
				);
			}

			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					grad_p_i -= gradRho;
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
					grad_p_i -= grad_p_j;
				);
			}		

			sum_grad_p_k += grad_p_i.squaredNorm();

			//////////////////////////////////////////////////////////////////////////
			// Compute factor as: factor_i = -1 / (a_ii * rho_i^2)
			// where a_ii is the diagonal entry of the linear system 
			// for the pressure A * p = source term
			//////////////////////////////////////////////////////////////////////////
			Real &factor = m_simulationData.getFactor(fluidModelIndex, i);
			if (sum_grad_p_k > m_eps)
				factor = static_cast<Real>(1.0) / (sum_grad_p_k);
			else
				factor = 0.0;
		}
	}
}

/** Compute rho_adv,i^(0) (see equation in Section 3.3 in [BK17])
  * using the velocities after the non-pressure forces were applied.
**/
void TimeStepDFSPH::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const Real h, const Real density0)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real &density = model->getDensity(i);
	Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	const Vector3r &xi = model->getPosition(i);
	const Vector3r &vi = model->getVelocity(i);
	Real delta = 0.0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors(
		const Vector3r & vj = fm_neighbor->getVelocity(neighborIndex);
		delta += (vi - vj).dot(sim->gradW(xi - xj));
		//delta += fm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
	);
	// assumes that all fluid particles have the same volume
	delta *= model->getVolume(i);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
			delta += bm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xi, vj);
			delta -= (vi - vj).dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xj, vj);
			delta += Vj * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}

	densityAdv = density / density0 + h*delta;
}

/** Compute rho_adv,i^(0) (see equation (9) in Section 3.2 [BK17])
  * using the velocities after the non-pressure forces were applied.
  */
void TimeStepDFSPH::computeDensityChange(const unsigned int fluidModelIndex, const unsigned int i, const Real h)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	const Vector3r &xi = model->getPosition(i);
	const Vector3r& vi = model->getVelocity(i);
	densityAdv = 0.0;
	unsigned int numNeighbors = 0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors(
		const Vector3r & vj = fm_neighbor->getVelocity(neighborIndex);
		densityAdv += (vi - vj).dot(sim->gradW(xi - xj));
	);
	// assumes that all fluid particles have the same volume
	densityAdv *= model->getVolume(i);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
			densityAdv += bm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xi, vj);
			densityAdv -= (vi - vj).dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xj, vj);
			densityAdv += Vj * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}
}

/** Compute pressure accelerations using the current pressure values of the particles
 */
void TimeStepDFSPH::computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, std::vector<std::vector<Real>>& pressure_rho2, const bool applyBoundaryForces)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
	ai.setZero();

	if (model->getParticleState(i) != ParticleState::Active)
		return;

	// p_rho2_i = (p_i / rho_i^2)
	const Real p_rho2_i = pressure_rho2[fluidModelIndex][i];
	const Vector3r &xi = model->getPosition(i);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors(			
		// p_rho2_j = (p_j / rho_j^2)
		const Real p_rho2_j = pressure_rho2[pid][neighborIndex];
		const Real pSum = p_rho2_i + fm_neighbor->getDensity0()/density0 * p_rho2_j;
		if (fabs(pSum) > m_eps)
		{
			const Vector3r grad_p_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
			ai += pSum * grad_p_j;		
		}
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (fabs(p_rho2_i) > m_eps)
	{
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			forall_boundary_neighbors(
				const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);

				const Vector3r a = (Real) 1.0 * p_rho2_i * grad_p_j;		
				ai += a;
				if (applyBoundaryForces)
					bm_neighbor->addForce(xj, -model->getMass(i) * a);
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		{
			forall_density_maps(
				const Vector3r a = (Real) 1.0 * p_rho2_i * gradRho;			
				ai += a;
				if (applyBoundaryForces)
					bm_neighbor->addForce(xj, -model->getMass(i) * a);
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		{
			forall_volume_maps(
				const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
				const Vector3r a = (Real) 1.0 * p_rho2_i * grad_p_j;		
				ai += a;

				if (applyBoundaryForces)
					bm_neighbor->addForce(xj, -model->getMass(i) * a);  
			);
		}
	}
}


Real TimeStepDFSPH::compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute A*p which is the change of the density when applying the 
	// pressure forces. 
	// \sum_j a_ij * p_j = h^2 \sum_j V_j (a_i - a_j) * gradW_ij
	// This is the RHS of Equation (12) in [BK17]
	//////////////////////////////////////////////////////////////////////////
	const Vector3r& xi = model->getPosition(i);
	const Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
	Real aij_pj = 0.0;

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors(
		const Vector3r & aj = m_simulationData.getPressureAccel(pid, neighborIndex);
		//aij_pj += fm_neighbor->getVolume(neighborIndex) * (ai - aj).dot(sim->gradW(xi - xj));
		aij_pj += (ai - aj).dot(sim->gradW(xi - xj));
	);
	// assumes that all fluid particles have the same volume
	aij_pj *= model->getVolume(i);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			aij_pj += bm_neighbor->getVolume(neighborIndex) * ai.dot(sim->gradW(xi - xj));
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			aij_pj -= ai.dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			aij_pj += Vj * ai.dot(sim->gradW(xi - xj));
		);
	}
	return aij_pj;
}


#endif
