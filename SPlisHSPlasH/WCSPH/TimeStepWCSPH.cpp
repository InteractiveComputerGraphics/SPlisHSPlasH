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
#include <SPlisHSPlasH/DynamicRigidBody.h>

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
	if (sim->getDynamicBoundarySimulator() != nullptr) {
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
	// Use dynamic boundary simulator and Akinci2012
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	// Compute density and v_s for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int r = 0; r < bm->numberOfParticles(); r++) {
				// compute v_s
				Vector3r gravForce = rb->getMass() * Vector3r(sim->getVecValue<Real>(Simulation::GRAVITATION));
				// the rb->getForce() contains only fluid-rigid forces
				Vector3r F_R = rb->getForce() + gravForce;
				bm->setV_s(r, rb->getVelocity() + dt * rb->getInvMass() * F_R + 
					(rb->getAngularVelocity() + rb->getInertiaTensorInverseW() * (dt * rb->getTorque() + dt * (rb->getInertiaTensorW() * rb->getAngularVelocity()).cross(rb->getAngularVelocity()))).cross(bm->getPosition(r)));

				// compute density
				Real particleDensity = 0;
				// iterate over all rigid bodies
				for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
					BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
					for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
						const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
						particleDensity += bm->getDensity0() * bm->getVolume(r) * sim->W(bm->getPosition0(r) - bm_neighbor->getPosition(k));
					}
				}
				bm->setDensity(r, particleDensity);
			}
		}
	}
	// compute source term s for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int r = 0; r < bm->numberOfParticles(); r++) {
				Real s = 0;
				// iterate over all rigid bodies except bm, since the divergence of particles in the same rigid body should be 0.
				for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
					if (boundaryPointSetIndex != pid) {
						BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
						for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
							const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
							// sum up divergence of v_s
							Real artificialVolumeK = (bm_neighbor->getDensity0() * bm_neighbor->getVolume(k)) / bm_neighbor->getDensity(k);
							s += artificialVolumeK * bm_neighbor->getDensity(k) * (bm_neighbor->getV_s(k) - bm->getV_s(r)).dot(sim->gradW(bm->getPosition(r) - bm_neighbor->getPosition(k)));
						}
					}
				}
				s += (bm->getDensity0() - bm->getDensity(r)) / dt;
				bm->setSourceTerm(r, s);
			}
		}
	}
}

void TimeStepWCSPH::performNeighborhoodSearch()
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

void TimeStepWCSPH::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepWCSPH::resize()
{
	m_simulationData.init();
}

