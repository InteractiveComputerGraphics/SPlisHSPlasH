#include "MicropolarModel_Bender2017.h"
#include <iostream>
#include "../TimeManager.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace GenParam;

int MicropolarModel_Bender2017::VISCOSITY_OMEGA = -1;
int MicropolarModel_Bender2017::INERTIA_INVERSE = -1;


MicropolarModel_Bender2017::MicropolarModel_Bender2017(FluidModel *model) :
	VorticityBase(model)
{
	m_omega.resize(model->numParticles(), Vector3r::Zero());
	m_angularAcceleration.resize(model->numParticles(), Vector3r::Zero());
	m_inertiaInverse = 0.5;
	m_viscosityOmega = 0.1;

	model->addField({ "angular velocity", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_omega[i][0]; }, true });
}

MicropolarModel_Bender2017::~MicropolarModel_Bender2017(void)
{
	m_model->removeFieldByName("angular velocity");

	m_omega.clear();
	m_angularAcceleration.clear();
}

void MicropolarModel_Bender2017::initParameters()
{
	VorticityBase::initParameters();

 	VISCOSITY_OMEGA = createNumericParameter("viscosityOmega", "Angular viscosity coefficient", &m_viscosityOmega);
 	setGroup(VISCOSITY_OMEGA, "Vorticity");
 	setDescription(VISCOSITY_OMEGA, "Viscosity coefficient for the angular velocity field.");
 	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_OMEGA));
 	rparam->setMinValue(0.0);


	INERTIA_INVERSE = createNumericParameter("inertiaInverse", "Inertia inverse", &m_inertiaInverse);
	setGroup(INERTIA_INVERSE, "Vorticity");
	setDescription(INERTIA_INVERSE, "Inverse microinertia used in the micropolar model.");
	rparam = static_cast<RealParameter*>(getParameter(INERTIA_INVERSE));
	rparam->setMinValue(0.0);
}


void MicropolarModel_Bender2017::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel *model = m_model;
	const Real density0 = model->getDensity0();

	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real invDt = static_cast<Real>(1.0) / dt;

	const Real nu_t = m_vorticityCoeff;
	const Real zeta = m_viscosityOmega;

	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;

	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 8.0;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			const Vector3r &omegai = m_omega[i];
			Vector3r &ai = m_model->getAcceleration(i);
			Vector3r &angAcceli = m_angularAcceleration[i];
			angAcceli.setZero();
			const Real density_i = m_model->getDensity(i);
			const Real density_i2 = density_i *density_i;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &vj = m_model->getVelocity(neighborIndex);
				const Vector3r &omegaj = m_omega[neighborIndex];

				// Viscosity
				const Real density_j = m_model->getDensity(neighborIndex);
				const Real density_j2 = density_j *density_j;
				const Vector3r xij = xi - xj;
				const Vector3r omegaij = omegai - omegaj;
				const Vector3r gradW = sim->gradW(xij);

 				// XSPH for angular velocity field
 				angAcceli -= invDt * m_inertiaInverse * zeta * (m_model->getMass(neighborIndex) / density_j) * omegaij * sim->W(xij);
				
				//// Viscosity
				//angAcceli += d * m_inertiaInverse * zeta * (m_model->getMass(neighborIndex) / density_i) * omegaij.dot(xij) / (xij.squaredNorm() + 0.01*h2) * gradW;

 				// difference curl 
 				ai += nu_t * 1.0/density_i * m_model->getMass(neighborIndex) * (omegaij.cross(gradW));
 				angAcceli += nu_t * 1.0/density_i * m_inertiaInverse * (m_model->getMass(neighborIndex) * (vi  - vj).cross(gradW));
			);

 			//////////////////////////////////////////////////////////////////////////
 			// Boundary
 			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
 					const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
 					const Vector3r &omegaj = Vector3r::Zero();//m_omega[neighborIndex];
 
 					// Viscosity
 					const Vector3r xij = xi - xj;
 					const Vector3r omegaij = omegai - omegaj;
 					const Vector3r gradW = sim->gradW(xij);
 
 					// XSPH for angular velocity field
 					//angAcceli -= invDt * m_inertiaInverse * zeta * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * omegaij * sim->W(xij);
 
 					// Viscosity
					//angAcceli += d * m_inertiaInverse * zeta * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * omegaij.dot(xij) / (xij.squaredNorm() + 0.01*h2) * gradW;
 
					// difference curl 
					const Vector3r a = nu_t * 1.0 / density_i * density0 * bm_neighbor->getVolume(neighborIndex) * (omegaij.cross(gradW));
					ai += a;
					angAcceli += nu_t * 1.0 / density_i * m_inertiaInverse * (density0 * bm_neighbor->getVolume(neighborIndex) * (vi - vj).cross(gradW));

					bm_neighbor->addForce(xj, -model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					Vector3r vj;
					bm_neighbor->getPointVelocity(xi, vj);
					const Vector3r omegaij = omegai;		// ToDo: omega_j

// 					// XSPH for angular velocity field
// 					angAcceli -= invDt * m_inertiaInverse * zeta * (density0 / density_i) * omegaij * rho;

					// difference curl 
					const Vector3r a = nu_t * density0 / density_i * (omegaij.cross(gradRho));
					ai += a;
					angAcceli += nu_t * density0 / density_i * m_inertiaInverse * ((vi - vj).cross(gradRho));

					bm_neighbor->addForce(xj, -model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					Vector3r vj;
					bm_neighbor->getPointVelocity(xj, vj);
					const Vector3r omegaij = omegai;		// ToDo: omega_j
					
					const Vector3r xij = xi - xj;
					const Vector3r gradW = sim->gradW(xij);

// 					// XSPH for angular velocity field
// 					angAcceli -= invDt * m_inertiaInverse * zeta * (density0 * Vj / density_i) * omegaij * sim->W(xij);

					// difference curl 
					const Vector3r a = nu_t * 1.0 / density_i * density0 * Vj * (omegaij.cross(gradW));
					ai += a;
					angAcceli += nu_t * 1.0 / density_i * m_inertiaInverse * (density0 * Vj * (vi - vj).cross(gradW));

					bm_neighbor->addForce(xj, -model->getMass(i) * a);
				);
			}
			angAcceli -= 2.0 * m_inertiaInverse * nu_t * omegai;
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_omega[i] += dt*m_angularAcceleration[i];
		}
	}
}


void MicropolarModel_Bender2017::reset()
{
	for (unsigned int i = 0; i < m_model->numActiveParticles(); i++)
		m_omega[i].setZero();
}

void SPH::MicropolarModel_Bender2017::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_omega[0]);
}

