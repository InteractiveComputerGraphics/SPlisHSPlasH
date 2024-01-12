#include "TimeStepIISPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataIISPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace std;

TimeStepIISPH::TimeStepIISPH() :
	TimeStep()
{
	m_simulationData.init();

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "a_ii", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getAii(fluidModelIndex, i); } });
		model->addField({ "d_ii", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDii(fluidModelIndex, i)[0]; } });
		model->addField({ "d_ij*p_j", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDij_pj(fluidModelIndex, i)[0]; } });
		model->addField({ "pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); }, true });
		//model->addField({ "last pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getLastPressure(fluidModelIndex, i); } });
		model->addField({ "advected density", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}
}

TimeStepIISPH::~TimeStepIISPH(void)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("a_ii");
		model->removeFieldByName("d_ii");
		model->removeFieldByName("d_ij*p_j");
		model->removeFieldByName("pressure");
		//model->removeFieldByName("last pressure");
		model->removeFieldByName("advected density");
		model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepIISPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		clearAccelerations(fluidModelIndex);

	sim->performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	precomputeValues();
#endif

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDensities(fluidModelIndex);

	sim->computeNonPressureForces();

	sim->updateTimeStepSize();

	// Solve density constraint	
	START_TIMING("predictAdvection");
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		predictAdvection(fluidModelIndex);
	STOP_TIMING_AVG;

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


void TimeStepIISPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
}

void TimeStepIISPH::predictAdvection(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = model->getDensity0();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Predict v_adv
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &vel = model->getVelocity(i);
			const Vector3r &accel = model->getAcceleration(i);
			Vector3r &dii = m_simulationData.getDii(fluidModelIndex, i);
			dii.setZero();

			if (model->getParticleState(i) == ParticleState::Active)
				vel += h * accel;

			// Compute d_ii
			const Real density = model->getDensity(i) / density0;
			const Real density2 = density*density;
			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				dii -= fm_neighbor->getVolume(neighborIndex) / density2 * sim->gradW(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					dii -= bm_neighbor->getVolume(neighborIndex) / density2 * sim->gradW(xi - xj);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					dii += 1.0 / density2 * gradRho;
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					dii -= Vj / density2 * sim->gradW(xi - xj);
				);
			}
		}
	}

	// Compute rho_adv
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Real density = model->getDensity(i) / density0;
			Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
			densityAdv = density;
			const Vector3r &xi = model->getPosition(i);
			const Vector3r &vi = model->getVelocity(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r &vj = fm_neighbor->getVelocity(neighborIndex);
				densityAdv += h*fm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
					densityAdv += h*bm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					Vector3r vj;
					bm_neighbor->getPointVelocity(xi, vj);
					densityAdv -= h*(vi - vj).dot(gradRho);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					Vector3r vj;
					bm_neighbor->getPointVelocity(xj, vj);
					densityAdv += h*Vj * (vi - vj).dot(sim->gradW(xi - xj));
				);
			}

			const Real &pressure = m_simulationData.getPressure(fluidModelIndex, i);
			Real &lastPressure = m_simulationData.getLastPressure(fluidModelIndex, i);
			lastPressure = static_cast<Real>(0.5)*pressure;

			// Compute a_ii
			Real &aii = m_simulationData.getAii(fluidModelIndex, i);
			aii = 0.0;

			const Vector3r &dii = m_simulationData.getDii(fluidModelIndex, i);

			const Real density2 = density * density;
			const Real dpi = model->getVolume(i) / density2;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				// Compute d_ji
				const Vector3r kernel = sim->gradW(xi - xj);
				const Vector3r dji = dpi * kernel;			
				aii += fm_neighbor->getVolume(neighborIndex) * (dii - dji).dot(kernel);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r kernel = sim->gradW(xi - xj);
					const Vector3r dji = dpi * kernel;			
					aii += bm_neighbor->getVolume(neighborIndex) * (dii - dji).dot(kernel);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					const Vector3r dji = -(1.0 / density2) * gradRho;
					aii -= (dii - dji).dot(gradRho);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r kernel = sim->gradW(xi - xj);
					const Vector3r dji = dpi * kernel;
					aii += Vj * (dii - dji).dot(kernel);
				);
			}
		}
	}
}

void TimeStepIISPH::pressureSolve()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

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
	INCREASE_COUNTER("IISPH - iterations", static_cast<Real>(m_iterations));
}

void TimeStepIISPH::pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	const Real density0 = model->getDensity0();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real omega = 0.5;

	// Compute dij_pj
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &dij_pj = m_simulationData.getDij_pj(fluidModelIndex, i);
			dij_pj.setZero();

			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Real densityj = fm_neighbor->getDensity(neighborIndex) / fm_neighbor->getDensity0();
				const Real densityj2 = densityj*densityj;

				dij_pj -= fm_neighbor->getVolume(neighborIndex) / densityj2 * m_simulationData.getLastPressure(pid, neighborIndex) * sim->gradW(xi - xj);
			);	
		}
	}

	// Compute new pressure
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real &pi = m_simulationData.getPressure(fluidModelIndex, i);
			pi = 0.0;

			const Real &aii = m_simulationData.getAii(fluidModelIndex, i);
			const Real density = model->getDensity(i) / density0;
			const Vector3r &xi = model->getPosition(i);

			const Real density2 = density*density;
			const Real dpi = model->getVolume(i) / density2;
			Real sum = 0.0;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r &d_jk_pk = m_simulationData.getDij_pj(pid, neighborIndex);

				// Compute \sum_{k \neq i} djk*pk
				// Compute d_ji
				const Vector3r kernel = sim->gradW(xi - xj);
				const Vector3r dji = dpi * kernel;
				const Vector3r d_ji_pi = dji * m_simulationData.getLastPressure(fluidModelIndex, i);

				// \sum ( mj * (\sum dij*pj - djj*pj - \sum_{k \neq i} djk*pk) * sim->gradW)
				sum += fm_neighbor->getVolume(neighborIndex) * (m_simulationData.getDij_pj(fluidModelIndex, i) - m_simulationData.getDii(pid, neighborIndex)*m_simulationData.getLastPressure(pid, neighborIndex) - (d_jk_pk - d_ji_pi)).dot(kernel);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					sum += bm_neighbor->getVolume(neighborIndex) * m_simulationData.getDij_pj(fluidModelIndex, i).dot(sim->gradW(xi - xj));
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					sum -= m_simulationData.getDij_pj(fluidModelIndex, i).dot(gradRho);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					sum += Vj * m_simulationData.getDij_pj(fluidModelIndex, i).dot(sim->gradW(xi - xj));
				);
			}

			const Real b = static_cast<Real>(1.0) - m_simulationData.getDensityAdv(fluidModelIndex, i);

			const Real &lastPi = m_simulationData.getLastPressure(fluidModelIndex, i);
			const Real denom = aii*h2;
			if (fabs(denom) > 1.0e-9)
				pi = max((static_cast<Real>(1.0) - omega)*lastPi + omega / denom * (b - h2*sum), static_cast<Real>(0.0));
			else
				pi = 0.0;

			if (pi != 0.0)
			{
				const Real newDensity = density0 * ((aii*pi + sum)*h2 - b) + density0;

				#pragma omp atomic
				avg_density_err += newDensity - density0;
			}
		}
	}

	for (int i = 0; i < (int)numParticles; i++)
	{
		const Real &pi = m_simulationData.getPressure(fluidModelIndex, i);
		Real &lastPi = m_simulationData.getLastPressure(fluidModelIndex, i);
		lastPi = pi;
	}

	avg_density_err /= numParticles;
}

void TimeStepIISPH::integration(const unsigned int fluidModelIndex)
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

void TimeStepIISPH::computePressureAccels(const unsigned int fluidModelIndex)
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
			const Real dpi = m_simulationData.getPressure(fluidModelIndex, i) / density2;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				// Pressure 
				const Real density_j = fm_neighbor->getDensity(neighborIndex) / fm_neighbor->getDensity0();
				const Real densityj2 = density_j*density_j;
				const Real dpj = m_simulationData.getPressure(pid, neighborIndex) / densityj2;
				ai -= fm_neighbor->getVolume(neighborIndex) * (dpi + fm_neighbor->getDensity0() / density0 * dpj) * sim->gradW(xi - xj);
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

void TimeStepIISPH::performNeighborhoodSearchSort()
{
	m_simulationData.performNeighborhoodSearchSort();
}

void TimeStepIISPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepIISPH::resize()
{
	m_simulationData.init();
}
