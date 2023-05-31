#include "TimeStepWCSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataWCSPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "Simulator/DynamicBoundarySimulator.h"
#include "SPlisHSPlasH/DynamicRigidBody.h"
#include "SPlisHSPlasH/StrongCouplingBoundarySolver.h"

using namespace SPH;
using namespace std;
using namespace GenParam;

int TimeStepWCSPH::STIFFNESS = -1;
int TimeStepWCSPH::EXPONENT = -1;


TimeStepWCSPH::TimeStepWCSPH() :
	TimeStep()
{
	m_simulationData.init();
	m_counter = 0;

	m_stiffness = 50.0;
	m_exponent = 7.0;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); } });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}
}

TimeStepWCSPH::~TimeStepWCSPH(void)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("pressure");
		model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepWCSPH::initParameters()
{
	TimeStep::initParameters();

	STIFFNESS = createNumericParameter("stiffness", "Stiffness", &m_stiffness);
	setGroup(STIFFNESS, "Simulation|WCSPH");
	setDescription(STIFFNESS, "Stiffness coefficient of EOS.");
	static_cast<RealParameter*>(getParameter(STIFFNESS))->setMinValue(static_cast<Real>(1e-6));

	EXPONENT = createNumericParameter("exponent", "Exponent (gamma)", &m_exponent);
	setGroup(EXPONENT, "Simulation|WCSPH");
	setDescription(EXPONENT, "Exponent of EOS.");
	static_cast<RealParameter*>(getParameter(EXPONENT))->setMinValue(static_cast<Real>(1e-6));

}

void TimeStepWCSPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	precomputeValues();
#endif

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	// Compute accelerations: a(t)
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		clearAccelerations(fluidModelIndex);
		computeDensities(fluidModelIndex);
	}
	sim->computeNonPressureForces();


	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const Real density0 = model->getDensity0();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				Real &density = model->getDensity(i);
				density = max(density, density0);
				m_simulationData.getPressure(fluidModelIndex, i) = m_stiffness * (pow(density / density0, m_exponent) - static_cast<Real>(1.0));
			}
		}

		computePressureAccels(fluidModelIndex);
	}

	// TODO: rigid-rigid contact forces
	if (sim->getDynamicBoundarySimulator() != nullptr && sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
		computeRigidRigidAccels();
	}

	sim->updateTimeStepSize();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static) 
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				if (model->getParticleState(i) == ParticleState::Active)
				{
					Vector3r &pos = model->getPosition(i);
					Vector3r &vel = model->getVelocity(i);
					Vector3r &accel = model->getAcceleration(i);
					accel += m_simulationData.getPressureAccel(fluidModelIndex, i);
					vel += accel * h;
					pos += vel * h;
				}
			}
		}
	}

	sim->emitParticles();
	sim->animateParticles();

	// Only for strong coupling method with BoundaryModel_Akinci2012
	if (sim->getDynamicBoundarySimulator() != nullptr && sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
		sim->getDynamicBoundarySimulator()->timeStepStrongCoupling();
	}

	// Compute new time	
	tm->setTime(tm->getTime() + h);
}


void TimeStepWCSPH::reset() {
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

void TimeStepWCSPH::computePressureAccels(const unsigned int fluidModelIndex) {
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	// Compute pressure forces
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++) {
			const Vector3r& xi = model->getPosition(i);
			const Real density_i = model->getDensity(i);

			Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
			ai.setZero();

			const Real dpi = m_simulationData.getPressure(fluidModelIndex, i) / (density_i * density_i);
			//////////////////////////////////////////////////////////////////////////
			// Fluid 
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				// Pressure 
				const Real density_j = fm_neighbor->getDensity(neighborIndex) * density0 / fm_neighbor->getDensity0();
			    const Real dpj = m_simulationData.getPressure(pid, neighborIndex) / (density_j * density_j);
			    ai -= density0 * fm_neighbor->getVolume(neighborIndex) * (dpi + dpj) * sim->gradW(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			const Real dpj = m_simulationData.getPressure(fluidModelIndex, i) / (density0 * density0);
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
				forall_boundary_neighbors(
					const Vector3r a = density0 * bm_neighbor->getVolume(neighborIndex) * (dpi + dpj) * sim->gradW(xi - xj);
				    ai -= a;
				    bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			} else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017) {
				forall_density_maps(
					const Vector3r a = -density0 * (dpi + dpj) * gradRho;
				    ai -= a;
				    bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			} else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019) {
				forall_volume_maps(
					const Vector3r a = density0 * Vj * (dpi + dpj) * sim->gradW(xi - xj);
				    ai -= a;
				    bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			}
		}
	}
}

void TimeStepWCSPH::computeRigidRigidAccels() {
	Simulation* sim = Simulation::getCurrent();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	StrongCouplingBoundarySolver* bs = StrongCouplingBoundarySolver::getCurrent();

	bs->computeDensityAndVolume();

	// compute density deviation
	int particleInContact = 0;
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		if (bm->getRigidBodyObject()->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static) reduction(+:particleInContact)
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						if (pid != boundaryPointSetIndex) {
							if (sim->numberOfNeighbors(boundaryPointSetIndex, pid, r) > 0) {
								particleInContact++;
								break;
							}
						}
					}
				}
			}
		}
	}
	if (particleInContact == 0) {
		return;
	}
	//std::cout << particleInContact << std::endl;
	Real avgDensityDeviation = 1.0;
	bs->computeDiagonalElement();
	bs->computeSourceTerm();
	int iterations = 0;
	// while density deviation too large
	while ((avgDensityDeviation >bs->getMaxDensityDeviation() && iterations < bs->getMaxIterations()) || iterations < bs->getMinIterations()) {
		//bs->computeSourceTermRHS();
		avgDensityDeviation = 0;
		// update pressure

		for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
			BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
			int bmIndex = boundaryPointSetIndex - nFluids;
			if (bm->getRigidBodyObject()->isDynamic()) {

				bs->computeSourceTermRHSForBody(boundaryPointSetIndex);

			    int numContactsBody = 0;
				#pragma omp parallel default(shared)
				{
                    #pragma omp for schedule(static) reduction(+:numContactsBody)
					for (int r = 0; r < bm->numberOfParticles(); r++) {
						for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
							if (pid != boundaryPointSetIndex) {
								if (sim->numberOfNeighbors(boundaryPointSetIndex, pid, r) > 0) {
									numContactsBody++;
									break;
								}
							}
						}
					}
				}
				if (numContactsBody == 0) {
					continue;
				}
				// beta_r_RJ
				Real relaxation = 0.5 / numContactsBody;

                #pragma omp parallel default(shared)
				{
                    #pragma omp for schedule(static) reduction(+:avgDensityDeviation)
					for (int r = 0; r < bm->numberOfParticles(); r++) {
						if (bs->getDensity(bmIndex, r) - bs->getRestDensity() > 1e-9 && std::abs(bs->getDiagonalElement(bmIndex, r)) > 1e-9) {
							Real pressureNextIter = bs->getPressure(bmIndex, r) + relaxation / bs->getDiagonalElement(bmIndex, r) * (bs->getSourceTerm(bmIndex, r) - bs->getSourceTermRHS(bmIndex, r));
							bs->setPressure(bmIndex, r, std::max(pressureNextIter, static_cast<Real>(0.0)));
							avgDensityDeviation -= (bs->getSourceTerm(bmIndex, r) - bs->getMinus_rho_div_v_s(bmIndex, r) - bs->getSourceTermRHS(bmIndex, r)) * dt;
						} else {
							bs->setPressure(bmIndex, r, static_cast<Real>(0.0));
						}
					}
				}
			}		
		}
		if (particleInContact == 0) {
			return;
		}
		avgDensityDeviation /= particleInContact;
		avgDensityDeviation /= bs->getRestDensity();
		avgDensityDeviation = std::abs(avgDensityDeviation);
		std::cout << avgDensityDeviation << std::endl;
		iterations++;
	}
	std::cout << "------------------------------" << std::endl;
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					const Vector3r f = -bs->getArtificialVolume(bmIndex, r) * bs->getPressureGrad(bmIndex, r);
					if (f != Vector3r::Zero()) {
						bm->addForce(bm->getPosition(r), f);					
					}
					bs->setPressure(bmIndex, r, static_cast<Real>(0.0));
				}
			}
		}
	}

	// solve the equation s = -rho * div¡¤v_rr w.r.t. pressure using relaxed jacobi

	//for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
	//	BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
	//	int bmIndex = boundaryPointSetIndex - nFluids;
	//	DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
	//	unsigned int numContacts = 0;
	//	if (bm->getRigidBodyObject()->isDynamic()) {
	//		bs->computeDensityAndVolume();
	//		// number of contacts
 //           #pragma omp parallel default(shared)
	//		{
 //               #pragma omp for schedule(static) reduction(+:numContacts)
	//			for (int r = 0; r < bm->numberOfParticles(); r++) {
	//				for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
	//					if (pid != boundaryPointSetIndex) {
	//						if (sim->numberOfNeighbors(boundaryPointSetIndex, pid, r) > 0) {
	//							numContacts++;
	//						}
	//					}
	//				}
	//			}
	//		}
	//		if(numContacts == 0){
	//			return;
	//		}
	//		// beta_r_RJ
	//		Real relaxation = 0.5/numContacts;
	//		//std::cout << numContacts << std::endl;
	//		bs->computeSourceTerm();
	//		bs->computeDiagonalElement();
	//		Real densityDeviationAvg = 1;
	//		unsigned int numParticles = 0;
	//		unsigned int i = 0;

	//		while (densityDeviationAvg > 0.002) {
	//			densityDeviationAvg = 0.0;
	//			numParticles = 0;
	//			bs->computeSourceTermRHS();
 //               #pragma omp parallel default(shared)
	//			{
 //                   #pragma omp for schedule(static) reduction(+:numParticles) reduction(+:densityDeviationAvg)
	//				for (int r = 0; r < bm->numberOfParticles(); r++) {
	//					if (bs->getDensity(bmIndex,r) - bs->getRestDensity(bmIndex) > 1e-9 && std::abs(bs->getDiagonalElement(bmIndex, r)) > 1e-9) {
	//						numParticles++;
	//						Real pressureNextIter = bs->getPressure(bmIndex, r) + relaxation / bs->getDiagonalElement(bmIndex, r) * (bs->getSourceTerm(bmIndex, r) - bs->getSourceTermRHS(bmIndex, r));
	//						//if (r == 105) {
	//						//	//if (bs->getSourceTerm(bmIndex, r) - bs->getSourceTermRHS(bmIndex, r) > 0) {
	//						//	//	std::cout << bs->getSourceTerm(bmIndex, r) - bs->getSourceTermRHS(bmIndex, r) << std::endl;
	//						//	//	std::cout << bs->getSourceTermRHS(bmIndex, r) << std::endl;
	//						//	//	std::cout << pressureNextIter << std::endl << std::endl;
	//						//	//}
	//						//	std::cout << bs->getSourceTerm(bmIndex, r) - bs->getSourceTermRHS(bmIndex, r) << std::endl;
	//						//	//std::cout << bs->getDiagonalElement(bmIndex, r) << std::endl;
	//						//	//std::cout << bs->getPressure(bmIndex, r) << std::endl;
	//						//	//pressureNextIter = 1.46289e+07;
	//						//}
	//						bs->setPressure(bmIndex, r, std::max(pressureNextIter, static_cast<Real>(0)));
	//						densityDeviationAvg -= (bs->getMinus_rho_div_v_s(bmIndex, r) + bs->getSourceTermRHS(bmIndex, r)) * dt;
	//					} else {
	//						bs->setPressure(bmIndex, r, 0);
	//					}
	//				}
	//			}
	//			if (numParticles == 0 || densityDeviationAvg < 0) {
	//				break;
	//			}
	//			densityDeviationAvg /= numParticles;
	//			densityDeviationAvg /= bs->getRestDensity(bmIndex);
	//			//std::cout << numParticles << std::endl;
	//			////densityErrorAvg /= numParticles;
	//			std::cout << densityDeviationAvg << std::endl;
	//			i++;
	//		}
	//		std::cout << "-------------------" << std::endl;
 //           #pragma omp parallel default(shared)
	//		{
 //               #pragma omp for schedule(static)  
	//			for (int r = 0; r < bm->numberOfParticles(); r++) {
	//				if (bs->getPressureGrad(bmIndex, r) != Vector3r::Zero()) {
	//					const Vector3r f = -bs->getArtificialVolume(bmIndex, r) * bs->getPressureGrad(bmIndex, r);
	//					if (bm->numberOfParticles() == 1) {
	//						bm->addForce(rb->getPosition(), f);
	//					} else {
	//						bm->addForce(bm->getPosition(r), f);
	//					}
	//					bs->setPressure(bmIndex, r, 0);
	//				} 
	//			}
	//		}
	//	}
	//}
}

void TimeStepWCSPH::performNeighborhoodSearch()
{
	if (Simulation::getCurrent()->zSortEnabled())
	{
		if (m_counter % 500 == 0)
		{
			Simulation::getCurrent()->performNeighborhoodSearchSort();
			StrongCouplingBoundarySolver::getCurrent()->performNeighborhoodSearchSort();
			m_simulationData.performNeighborhoodSearchSort();
		}
		m_counter++;
	}

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepWCSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepWCSPH::resize()
{
	m_simulationData.init();
}

