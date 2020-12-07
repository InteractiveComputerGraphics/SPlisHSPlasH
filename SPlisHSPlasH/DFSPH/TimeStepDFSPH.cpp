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
	m_counter = 0;
	m_iterationsV = 0;
	m_enableDivergenceSolver = true;
	m_maxIterationsV = 100;
	m_maxErrorV = static_cast<Real>(0.1);

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "factor", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getFactor(fluidModelIndex, i); } });
		model->addField({ "advected density", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "kappa", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getKappa(fluidModelIndex, i); }, true });
		model->addField({ "kappa_v", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getKappaV(fluidModelIndex, i); }, true });
	}
}

TimeStepDFSPH::~TimeStepDFSPH(void)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("factor");
		model->removeFieldByName("advected density");
		model->removeFieldByName("kappa");
		model->removeFieldByName("kappa_v");
	}
}

void TimeStepDFSPH::initParameters()
{
	TimeStep::initParameters();

	SOLVER_ITERATIONS_V = createNumericParameter("iterationsV", "Iterations (divergence)", &m_iterationsV);
	setGroup(SOLVER_ITERATIONS_V, "DFSPH");
	setDescription(SOLVER_ITERATIONS_V, "Iterations required by the divergence solver.");
	getParameter(SOLVER_ITERATIONS_V)->setReadOnly(true);

	MAX_ITERATIONS_V = createNumericParameter("maxIterationsV", "Max. iterations (divergence)", &m_maxIterationsV);
	setGroup(MAX_ITERATIONS_V, "DFSPH");
	setDescription(MAX_ITERATIONS_V, "Maximal number of iterations of the divergence solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS_V))->setMinValue(1);

	MAX_ERROR_V = createNumericParameter("maxErrorV", "Max. divergence error(%)", &m_maxErrorV);
	setGroup(MAX_ERROR_V, "DFSPH");
	setDescription(MAX_ERROR_V, "Maximal divergence error (%).");
	static_cast<RealParameter*>(getParameter(MAX_ERROR_V))->setMinValue(static_cast<Real>(1e-6));

	USE_DIVERGENCE_SOLVER = createBoolParameter("enableDivergenceSolver", "Enable divergence solver", &m_enableDivergenceSolver);
	setGroup(USE_DIVERGENCE_SOLVER, "DFSPH");
	setDescription(USE_DIVERGENCE_SOLVER, "Turn divergence solver on/off.");
}

void TimeStepDFSPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();
	const unsigned int nModels = sim->numberOfFluidModels();

	performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	precomputeValues();
#endif

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDensities(fluidModelIndex);

	START_TIMING("computeDFSPHFactor");
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDFSPHFactor(fluidModelIndex);
	STOP_TIMING_AVG;

	if (m_enableDivergenceSolver)
	{
		START_TIMING("divergenceSolve");
		divergenceSolve();
		STOP_TIMING_AVG
	}
	else
		m_iterationsV = 0;

	// Compute accelerations: a(t)
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		clearAccelerations(fluidModelIndex);

	sim->computeNonPressureForces();

	sim->updateTimeStepSize();

	// compute new velocities only considering non-pressure forces
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

	// compute final positions
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

	sim->emitParticles();
	sim->animateParticles();

	// Compute new time	
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

#ifdef USE_WARMSTART	
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
		warmstartPressureSolve(fluidModelIndex);
#endif

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
				computeDensityAdv(fluidModelIndex, i, numParticles, h, density0);
				m_simulationData.getFactor(fluidModelIndex, i) *= invH2;
//#ifdef USE_WARMSTART
//				m_simulationData.getKappa(fluidModelIndex, i) = 0.0;
//#endif
			}
		}
	}

	m_iterations = 0;

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////
	
	Real avg_density_err = 0.0;
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

	INCREASE_COUNTER("DFSPH - iterations", static_cast<Real>(m_iterations));

#ifdef USE_WARMSTART
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();

		//////////////////////////////////////////////////////////////////////////
		// Multiply by h^2, the time step size has to be removed 
		// to make the stiffness value independent 
		// of the time step size
		//////////////////////////////////////////////////////////////////////////
		for (int i = 0; i < numParticles; i++)
			m_simulationData.getKappa(fluidModelIndex, i) *= h2;
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


#ifdef USE_WARMSTART_V
	for(unsigned int fluidModelIndex =0; fluidModelIndex < nFluids; fluidModelIndex++)
		warmstartDivergenceSolve(fluidModelIndex);
#endif

	//////////////////////////////////////////////////////////////////////////
	// Compute velocity of density change
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
				computeDensityChange(fluidModelIndex, i, h);
				m_simulationData.getFactor(fluidModelIndex, i) *= invH;

//#ifdef USE_WARMSTART_V
//				m_simulationData.getKappaV(fluidModelIndex, i) = 0.0;
//#endif
			}
		}
	}

	m_iterationsV = 0;

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////
	
	Real avg_density_err = 0.0;
	bool chk = false;

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

	//////////////////////////////////////////////////////////////////////////
	// Multiply by h, the time step size has to be removed 
	// to make the stiffness value independent 
	// of the time step size
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();

#ifdef USE_WARMSTART_V
		for (int i = 0; i < numParticles; i++)
			m_simulationData.getKappaV(fluidModelIndex, i) *= h;
#endif

		for (int i = 0; i < numParticles; i++)
		{
			m_simulationData.getFactor(fluidModelIndex, i) *= h;
		}
	}
}


#ifdef USE_AVX

 void TimeStepDFSPH::computeDFSPHFactor(const unsigned int fluidModelIndex)
 {
 	//////////////////////////////////////////////////////////////////////////
 	// Init parameters
 	//////////////////////////////////////////////////////////////////////////
 
 	Simulation *sim = Simulation::getCurrent();
 	const unsigned int nFluids = sim->numberOfFluidModels();
 	FluidModel *model = sim->getFluidModel(fluidModelIndex);
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
 			// Compute gradient dp_i/dx_j * (1/k)  and dp_j/dx_j * (1/k)
 			//////////////////////////////////////////////////////////////////////////
 			const Vector3r &xi = model->getPosition(i);
 
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
				const Vector3f8 &gradC_j = V_gradW;
 				sum_grad_p_k_avx += gradC_j.squaredNorm();
 				grad_p_i_avx += gradC_j;
 			);
 
 			//////////////////////////////////////////////////////////////////////////
 			// Boundary
 			//////////////////////////////////////////////////////////////////////////
 			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
 			{
 				forall_boundary_neighbors_avx(
 					const Scalarf8 V_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
 					const Vector3f8 gradC_j = CubicKernel_AVX::gradW(xj_avx - xi_avx) * V_avx;
 					grad_p_i_avx -= gradC_j;
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
 			// Compute pressure stiffness denominator
 			//////////////////////////////////////////////////////////////////////////
 			Real &factor = m_simulationData.getFactor(fluidModelIndex, i);
 			if (sum_grad_p_k > m_eps)
 				factor = -static_cast<Real>(1.0) / (sum_grad_p_k);
 			else
 				factor = 0.0;
 		}
 	}
 }

#ifdef USE_WARMSTART
void TimeStepDFSPH::warmstartPressureSolve(const unsigned int fluidModelIndex)
{
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h * h;
	const Real invH = static_cast<Real>(1.0) / h;
	const Real invH2 = static_cast<Real>(1.0) / h2;
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	const Scalarf8 h_avx(h);
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Scalarf8 invH_avx(invH);

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Divide by h^2, the time step size has been removed in 
		// the last step to make the stiffness value independent 
		// of the time step size
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			//m_simulationData.getKappa(fluidModelIndex, i) = max(m_simulationData.getKappa(fluidModelIndex, i)*invH2, -static_cast<Real>(0.5) * density0*density0);
			computeDensityAdv(fluidModelIndex, i, numParticles, h, density0);
			if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 1.0)
				m_simulationData.getKappa(fluidModelIndex, i) = static_cast<Real>(0.5) * max(m_simulationData.getKappa(fluidModelIndex, i), static_cast<Real>(-0.00025)) * invH2;
			else
				m_simulationData.getKappa(fluidModelIndex, i) = 0.0;
		}

		//////////////////////////////////////////////////////////////////////////
		// Predict v_adv with external velocities
		////////////////////////////////////////////////////////////////////////// 

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
			{
				m_simulationData.getKappa(fluidModelIndex, i) = 0.0;
				continue;
			}

			//if (m_simulationData.getDensityAdv(fluidModelIndex, i) > density0)
			{
				const Real ki = m_simulationData.getKappa(fluidModelIndex, i);
				const Vector3r &xi = model->getPosition(i);
				Vector3r &vi = model->getVelocity(i);

				Scalarf8 ki_avx(ki);
				Vector3f8 xi_avx(xi);
				Vector3f8 delta_vi;
				delta_vi.setZero();

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				forall_fluid_neighbors_avx_nox(
					compute_xj(fm_neighbor, pid);
					compute_Vj(fm_neighbor);
					compute_Vj_gradW();
					const Scalarf8 densityFrac_avx(fm_neighbor->getDensity0() / density0);
					const Scalarf8 kj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &m_simulationData.getKappa(pid, 0), count);
					const Scalarf8 kSum_avx = ki_avx + densityFrac_avx * kj_avx;

					// Directly update velocities instead of storing pressure accelerations
					delta_vi += V_gradW * (h_avx * kSum_avx);			// ki, kj already contain inverse density	
				);

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (fabs(ki) > m_eps)
				{
					if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
					{
						const Scalarf8 mi_avx(model->getMass(i));
						forall_boundary_neighbors_avx(
							const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);

							// Directly update velocities instead of storing pressure accelerations
							const Vector3f8 velChange = -CubicKernel_AVX::gradW(xj_avx - xi_avx) * (Vj_avx * Scalarf8(1.0f*h) * ki_avx);			// ki, kj already contain inverse density	
							delta_vi += velChange;

							bm_neighbor->addForce(xj_avx, -velChange * (mi_avx*invH_avx), count);
						);
					}
				}
				vi[0] += delta_vi.x().reduce();
				vi[1] += delta_vi.y().reduce();
				vi[2] += delta_vi.z().reduce();

				if (fabs(ki) > m_eps)
				{
					if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
					{
						forall_density_maps(
							const Vector3r velChange = -h * ki * gradRho;				// kj already contains inverse density
						vi += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
					else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
					{
						forall_volume_maps(
							const Vector3r velChange = h * ki * Vj * sim->gradW(xi - xj);				// kj already contains inverse density
							vi += velChange;
							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
				}
			}
		}
	}
}
#endif

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
	const Scalarf8 invH_avx(invH);
	const Scalarf8 h_avx(h);

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure forces
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
				continue;

			//////////////////////////////////////////////////////////////////////////
			// Evaluate rhs
			//////////////////////////////////////////////////////////////////////////
			const Real b_i = m_simulationData.getDensityAdv(fluidModelIndex, i) - static_cast<Real>(1.0);
			const Real ki = b_i*m_simulationData.getFactor(fluidModelIndex, i);
#ifdef USE_WARMSTART
			m_simulationData.getKappa(fluidModelIndex, i) += ki;
#endif

			Vector3r &vi = model->getVelocity(i);
			const Vector3r &xi = model->getPosition(i);

			Scalarf8 ki_avx(ki);
			Vector3f8 xi_avx(xi);
			Vector3f8 delta_vi;
			delta_vi.setZero();


			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_avx_nox(
				compute_xj(fm_neighbor, pid);
				compute_Vj(fm_neighbor);
				compute_Vj_gradW();
				const Scalarf8 densityFrac_avx(fm_neighbor->getDensity0() / density0);
				const Scalarf8 densityAdvj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &m_simulationData.getDensityAdv(pid, 0), count);
				const Scalarf8 factorj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &m_simulationData.getFactor(pid, 0), count);
					
				const Scalarf8 b_j_avx = densityAdvj_avx - Scalarf8(1.0f);
				const Scalarf8 kj_avx = b_j_avx * factorj_avx;
				const Scalarf8 kSum_avx = ki_avx + densityFrac_avx * kj_avx;

				// Directly update velocities instead of storing pressure accelerations
				delta_vi += V_gradW * (h_avx * kSum_avx);			// ki, kj already contain inverse density	
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (fabs(ki) > m_eps)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					const Scalarf8 mi_avx(model->getMass(i));
					forall_boundary_neighbors_avx(
						const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);

						// Directly update velocities instead of storing pressure accelerations
						const Vector3f8 velChange = -CubicKernel_AVX::gradW(xj_avx - xi_avx) * (Vj_avx * Scalarf8(1.0f*h) * ki_avx);			// ki, kj already contain inverse density	
						delta_vi += velChange;

						bm_neighbor->addForce(xj_avx, -velChange * (mi_avx*invH_avx), count);
					);
				}
			}	
			
			vi[0] += delta_vi.x().reduce();
			vi[1] += delta_vi.y().reduce();
			vi[2] += delta_vi.z().reduce();
			
			if (fabs(ki) > m_eps)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						const Vector3r velChange = -h * ki * gradRho;				// kj already contains inverse density
						vi += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						const Vector3r velChange = h * ki * Vj * sim->gradW(xi - xj);				// kj already contains inverse density
						vi += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
			}
		}
	
		//////////////////////////////////////////////////////////////////////////
		// Update rho_adv and density error
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			computeDensityAdv(fluidModelIndex, i, numParticles, h, density0);

			density_error += m_simulationData.getDensityAdv(fluidModelIndex, i) - static_cast<Real>(1.0);
		}
	}
	avg_density_err = density0 * density_error / numParticles;
}

#ifdef USE_WARMSTART_V
void TimeStepDFSPH::warmstartDivergenceSolve(const unsigned int fluidModelIndex)
{
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Scalarf8 invH_avx(invH);
	const Scalarf8 h_avx(h);

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Divide by h^2, the time step size has been removed in 
		// the last step to make the stiffness value independent 
		// of the time step size
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			computeDensityChange(fluidModelIndex, i, h);
			if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 0.0)
				m_simulationData.getKappaV(fluidModelIndex, i) = static_cast<Real>(0.5) * max(m_simulationData.getKappaV(fluidModelIndex, i), static_cast<Real>(-0.5)) * invH;
			else
				m_simulationData.getKappaV(fluidModelIndex, i) = 0.0;
		}

		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
			{
				m_simulationData.getKappaV(fluidModelIndex, i) = 0.0;
				continue;
			}

			//if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 0.0)
			{
				const Real ki = m_simulationData.getKappaV(fluidModelIndex, i);
				const Vector3r &xi = model->getPosition(i);
				Vector3r &vi = model->getVelocity(i);

				Scalarf8 ki_avx(ki);
				Vector3f8 xi_avx(xi);
				Vector3f8 delta_vi;
				delta_vi.setZero();

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				forall_fluid_neighbors_avx_nox(
					compute_xj(fm_neighbor, pid);
					compute_Vj(fm_neighbor);
					compute_Vj_gradW();
					const Scalarf8 densityFrac_avx(fm_neighbor->getDensity0() / density0);
					const Scalarf8 kj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &m_simulationData.getKappaV(pid, 0), count);;
					const Scalarf8 kSum_avx = ki_avx + densityFrac_avx * kj_avx;

					// Directly update velocities instead of storing pressure accelerations
					delta_vi += V_gradW * ( h_avx * kSum_avx);			// ki, kj already contain inverse density	
				);

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (fabs(ki) > m_eps)
				{
					if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
					{
						const Scalarf8 mi_avx(model->getMass(i));
						forall_boundary_neighbors_avx(
							const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);

							// Directly update velocities instead of storing pressure accelerations
							const Vector3f8 velChange = -CubicKernel_AVX::gradW(xj_avx - xi_avx) * (Vj_avx * h_avx * ki_avx);			// ki, kj already contain inverse density	
							delta_vi += velChange;

							bm_neighbor->addForce(xj_avx, -velChange * (mi_avx*invH_avx), count);
						);
					}
				}
				vi[0] += delta_vi.x().reduce();
				vi[1] += delta_vi.y().reduce();
				vi[2] += delta_vi.z().reduce();

				if (fabs(ki) > m_eps)
				{
					if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
					{
						forall_density_maps(
							const Vector3r velChange = -h * ki * gradRho;				// kj already contains inverse density
							vi += velChange;
							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
					else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
					{
						forall_volume_maps(
							const Vector3r velChange = h * ki * Vj * sim->gradW(xi - xj);				// kj already contains inverse density
							vi += velChange;
							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
				}
			}
		}
	}
}
#endif

void TimeStepDFSPH::divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	Real density_error = 0.0;
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Scalarf8 invH_avx(invH);
	const Scalarf8 h_avx(h);

	//////////////////////////////////////////////////////////////////////////
	// Perform Jacobi iteration over all blocks
	//////////////////////////////////////////////////////////////////////////	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
				continue;

			//////////////////////////////////////////////////////////////////////////
			// Evaluate rhs
			//////////////////////////////////////////////////////////////////////////
			const Real b_i = m_simulationData.getDensityAdv(fluidModelIndex, i);
			const Real ki = b_i*m_simulationData.getFactor(fluidModelIndex, i);
#ifdef USE_WARMSTART_V
			m_simulationData.getKappaV(fluidModelIndex, i) += ki;
#endif

			Vector3r &vi = model->getVelocity(i);
			const Vector3r &xi = model->getPosition(i);

			Scalarf8 ki_avx(ki);
			Vector3f8 xi_avx(xi);
			Vector3f8 delta_vi;
			delta_vi.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////			
			forall_fluid_neighbors_avx_nox(
				compute_xj(fm_neighbor, pid);
				compute_Vj(fm_neighbor);
				compute_Vj_gradW();
				const Scalarf8 densityFrac_avx(fm_neighbor->getDensity0() / density0);
				const Scalarf8 densityAdvj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &m_simulationData.getDensityAdv(pid, 0), count);
				const Scalarf8 factorj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &m_simulationData.getFactor(pid, 0), count);

				const Scalarf8 kj_avx = densityAdvj_avx * factorj_avx;
				const Scalarf8 kSum_avx = ki_avx + densityFrac_avx * kj_avx;

				// Directly update velocities instead of storing pressure accelerations
				delta_vi += V_gradW * (h_avx * kSum_avx);			// ki, kj already contain inverse density	
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (fabs(ki) > m_eps)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					const Scalarf8 mi_avx(model->getMass(i));
					forall_boundary_neighbors_avx(
						const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);

						// Directly update velocities instead of storing pressure accelerations
						const Vector3f8 velChange = -CubicKernel_AVX::gradW(xj_avx - xi_avx) * (Vj_avx * h_avx * ki_avx);			// ki, kj already contain inverse density	
						delta_vi += velChange;

						bm_neighbor->addForce(xj_avx, -velChange * (mi_avx*invH_avx), count);
					);
				}
			}

			vi[0] += delta_vi.x().reduce();
			vi[1] += delta_vi.y().reduce();
			vi[2] += delta_vi.z().reduce();

			if (fabs(ki) > m_eps)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						const Vector3r velChange = -h * ki * gradRho;				// kj already contains inverse density
						vi += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						const Vector3r velChange = h * ki * Vj * sim->gradW(xi - xj);				// kj already contains inverse density
						vi += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
			}
		}
		//////////////////////////////////////////////////////////////////////////
		// Update rho_adv and density error
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			computeDensityChange(fluidModelIndex, i, h);
			density_error += m_simulationData.getDensityAdv(fluidModelIndex, i);
		}
	}
	avg_density_err = density0 * density_error/numParticles;
}

void TimeStepDFSPH::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const int numParticles, const Real h, const Real density0)
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
	densityAdv = max(densityAdv, static_cast<Real>(1.0));
}

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

	// only correct positive divergence
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

	densityAdv = max(densityAdv, static_cast<Real>(0.0));

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
			// Compute gradient dp_i/dx_j * (1/k)  and dp_j/dx_j * (1/k)
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
			// Compute pressure stiffness denominator
			//////////////////////////////////////////////////////////////////////////
			Real &factor = m_simulationData.getFactor(fluidModelIndex, i);
			if (sum_grad_p_k > m_eps)
				factor = -static_cast<Real>(1.0) / (sum_grad_p_k);
			else
				factor = 0.0;
		}
	}
}

#ifdef USE_WARMSTART
void TimeStepDFSPH::warmstartPressureSolve(const unsigned int fluidModelIndex)
{
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real invH = static_cast<Real>(1.0) / h;
	const Real invH2 = static_cast<Real>(1.0) / h2;
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Divide by h^2, the time step size has been removed in 
		// the last step to make the stiffness value independent 
		// of the time step size
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			//m_simulationData.getKappa(fluidModelIndex, i) = max(m_simulationData.getKappa(fluidModelIndex, i)*invH2, -static_cast<Real>(0.5) * density0*density0);
			computeDensityAdv(fluidModelIndex, i, numParticles, h, density0);
			if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 1.0)
				m_simulationData.getKappa(fluidModelIndex, i) = static_cast<Real>(0.5) * max(m_simulationData.getKappa(fluidModelIndex, i), static_cast<Real>(-0.00025)) * invH2;
			else
				m_simulationData.getKappa(fluidModelIndex, i) = 0.0;
		}

		//////////////////////////////////////////////////////////////////////////
		// Predict v_adv with external velocities
		////////////////////////////////////////////////////////////////////////// 

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
			{
				m_simulationData.getKappa(fluidModelIndex, i) = 0.0;
				continue;
			}

			//if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 1.0)
			{
				Vector3r &vel = model->getVelocity(i);
				const Real ki = m_simulationData.getKappa(fluidModelIndex, i);
				const Vector3r &xi = model->getPosition(i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				forall_fluid_neighbors(
					const Real kj = m_simulationData.getKappa(pid, neighborIndex);

					const Real kSum = (ki + fm_neighbor->getDensity0() / density0 * kj);
					if (fabs(kSum) > m_eps)
					{
						const Vector3r grad_p_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
						vel -= h * kSum * grad_p_j;					// ki, kj already contain inverse density
					}
				)

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (fabs(ki) > m_eps)
				{
					if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
					{
						forall_boundary_neighbors(
							const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
							const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
							vel += velChange;
							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
					else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
					{
						forall_density_maps(
							const Vector3r velChange = -h * (Real) 1.0 * ki * gradRho;				// kj already contains inverse density
							vel += velChange;
							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
					else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
					{
						forall_volume_maps(
							const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
							const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
							vel += velChange;

							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
				}
			}
		}
	}
}
#endif 

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
		// Compute pressure forces
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
				continue;

			//////////////////////////////////////////////////////////////////////////
			// Evaluate rhs
			//////////////////////////////////////////////////////////////////////////
			const Real b_i = m_simulationData.getDensityAdv(fluidModelIndex, i) - static_cast<Real>(1.0);
			const Real ki = b_i*m_simulationData.getFactor(fluidModelIndex, i);
#ifdef USE_WARMSTART
			m_simulationData.getKappa(fluidModelIndex, i) += ki;
#endif

			Vector3r &v_i = model->getVelocity(i);
			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(				
				const Real b_j = m_simulationData.getDensityAdv(pid, neighborIndex) - static_cast<Real>(1.0);
				const Real kj = b_j*m_simulationData.getFactor(pid, neighborIndex);
				const Real kSum = ki + fm_neighbor->getDensity0()/density0 * kj;
				if (fabs(kSum) > m_eps)
				{
					const Vector3r grad_p_j = -fm_neighbor->getVolume(neighborIndex) *sim->gradW(xi - xj);

					// Directly update velocities instead of storing pressure accelerations
					v_i -= h * kSum * grad_p_j;			// ki, kj already contain inverse density						
				}
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (fabs(ki) > m_eps)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);

						// Directly update velocities instead of storing pressure accelerations
						const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
						v_i += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						const Vector3r velChange = -h * (Real) 1.0 * ki * gradRho;				// kj already contains inverse density
						v_i += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
						const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
						v_i += velChange;

						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
			}
		}


		//////////////////////////////////////////////////////////////////////////
		// Update rho_adv and density error
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			computeDensityAdv(fluidModelIndex, i, numParticles, h, density0);

			density_error += density0 * m_simulationData.getDensityAdv(fluidModelIndex, i) - density0;
		}
	}

	avg_density_err = density_error / numParticles;
}

#ifdef USE_WARMSTART_V
void TimeStepDFSPH::warmstartDivergenceSolve(const unsigned int fluidModelIndex)
{
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();


	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Divide by h^2, the time step size has been removed in 
		// the last step to make the stiffness value independent 
		// of the time step size
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			//m_simulationData.getKappaV(fluidModelIndex, i) = static_cast<Real>(0.5)*max(m_simulationData.getKappaV(fluidModelIndex, i)*invH, -static_cast<Real>(0.5) * density0*density0);
			computeDensityChange(fluidModelIndex, i, h);
			if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 0.0)				
				m_simulationData.getKappaV(fluidModelIndex, i) = static_cast<Real>(0.5) * max(m_simulationData.getKappaV(fluidModelIndex, i), static_cast<Real>(-0.5)) * invH;
			else
				m_simulationData.getKappaV(fluidModelIndex, i) = 0.0;
		}

		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
			{
				m_simulationData.getKappaV(fluidModelIndex, i) = 0.0;
				continue;
			}

			//if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 0.0)
			{
				Vector3r &vel = model->getVelocity(i);
				const Real ki = m_simulationData.getKappaV(fluidModelIndex, i);
				const Vector3r &xi = model->getPosition(i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				forall_fluid_neighbors(
					const Real kj = m_simulationData.getKappaV(pid, neighborIndex);

					const Real kSum = (ki + fm_neighbor->getDensity0() / density0 * kj);
					if (fabs(kSum) > m_eps)
					{
						const Vector3r grad_p_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
						vel -= h * kSum * grad_p_j;					// ki, kj already contain inverse density
					}
				)

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (fabs(ki) > m_eps)
				{
					if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
					{
						forall_boundary_neighbors(
							const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
							const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
							vel += velChange;
							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
					else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
					{
						forall_density_maps(
							const Vector3r velChange = -h * (Real) 1.0 * ki * gradRho;				// kj already contains inverse density
							vel += velChange;
							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
					else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
					{
						forall_volume_maps(
							const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
							const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
							vel += velChange;
	
							bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
						);
					}
				}
			}
		}
	}
}
#endif

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

	//////////////////////////////////////////////////////////////////////////
	// Perform Jacobi iteration over all blocks
	//////////////////////////////////////////////////////////////////////////	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
				continue;

			//////////////////////////////////////////////////////////////////////////
			// Evaluate rhs
			//////////////////////////////////////////////////////////////////////////
			const Real b_i = m_simulationData.getDensityAdv(fluidModelIndex, i);
			const Real ki = b_i*m_simulationData.getFactor(fluidModelIndex, i);
#ifdef USE_WARMSTART_V
			m_simulationData.getKappaV(fluidModelIndex, i) += ki;
#endif

			Vector3r &v_i = model->getVelocity(i);

			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Real b_j = m_simulationData.getDensityAdv(pid, neighborIndex);
				const Real kj = b_j*m_simulationData.getFactor(pid, neighborIndex);

				const Real kSum = ki + fm_neighbor->getDensity0() / density0 * kj;
				if (fabs(kSum) > m_eps)
				{
					const Vector3r grad_p_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
					v_i -= h * kSum * grad_p_j;			// ki, kj already contain inverse density
				}
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (fabs(ki) > m_eps)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
						const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
						v_i += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						const Vector3r velChange = -h * (Real) 1.0 * ki * gradRho;				// kj already contains inverse density
						v_i += velChange;
						bm_neighbor->addForce(xj, -model->getMass(i) * velChange * invH);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
						const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
						v_i += velChange;

						bm_neighbor->addForce(xj, - model->getMass(i) * velChange * invH);
					);
				}
			}
		}

		//////////////////////////////////////////////////////////////////////////
		// Update rho_adv and density error
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			computeDensityChange(fluidModelIndex, i, h);
			density_error += density0 * m_simulationData.getDensityAdv(fluidModelIndex, i);
		}
	}
	avg_density_err = density_error/numParticles;
}

void TimeStepDFSPH::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const int numParticles, const Real h, const Real density0)
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
		const Vector3r &vj = fm_neighbor->getVelocity(neighborIndex);
		delta += fm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
	)

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
	densityAdv = max(densityAdv, static_cast<Real>(1.0));
}

void TimeStepDFSPH::computeDensityChange(const unsigned int fluidModelIndex, const unsigned int i, const Real h)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	const Vector3r &xi = model->getPosition(i);
	const Vector3r &vi = model->getVelocity(i);
	densityAdv = 0.0;
	unsigned int numNeighbors = 0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors(
		const Vector3r &vj = fm_neighbor->getVelocity(neighborIndex);
		densityAdv += fm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
	);

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

	// only correct positive divergence
	densityAdv = max(densityAdv, static_cast<Real>(0.0));

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
}

#endif

void TimeStepDFSPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
	m_iterations = 0;
	m_iterationsV = 0;
}

void TimeStepDFSPH::performNeighborhoodSearch()
{
	if (Simulation::getCurrent()->zSortEnabled())
	{
		if (m_counter % 500 == 0)
		{
			Simulation::getCurrent()->performNeighborhoodSearchSort();
			m_simulationData.performNeighborhoodSearchSort();
		}
		m_counter++;
	}

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepDFSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepDFSPH::resize()
{
	m_simulationData.init();
}
