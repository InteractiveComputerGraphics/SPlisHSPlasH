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
	m_inertiaInverse = static_cast<Real>(0.5);
	m_viscosityOmega = static_cast<Real>(0.1);

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
 	setGroup(VISCOSITY_OMEGA, "Fluid Model|Vorticity");
 	setDescription(VISCOSITY_OMEGA, "Viscosity coefficient for the angular velocity field.");
 	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_OMEGA));
 	rparam->setMinValue(0.0);


	INERTIA_INVERSE = createNumericParameter("inertiaInverse", "Inertia inverse", &m_inertiaInverse);
	setGroup(INERTIA_INVERSE, "Fluid Model|Vorticity");
	setDescription(INERTIA_INVERSE, "Inverse microinertia used in the micropolar model.");
	rparam = static_cast<RealParameter*>(getParameter(INERTIA_INVERSE));
	rparam->setMinValue(0.0);
}

#ifdef USE_AVX

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
	const Scalarf8 density0_avx(density0);

	const Scalarf8 factor_avx(invDt * m_inertiaInverse * zeta *density0);

	//Real d = 10.0;
	//if (sim->is2DSimulation())
	//	d = 8.0;

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

			const Vector3f8 xi_avx(xi);
			const Vector3f8 vi_avx(vi);
			const Scalarf8 mi_avx(m_model->getMass(i));
			const Vector3f8 omegai_avx(omegai);
			const Scalarf8 density_i_avx(density_i);
			const Scalarf8 nut_density_i(nu_t / density_i);
			const Scalarf8 nut_density_i_intertiaInverse(nu_t / density_i * m_inertiaInverse);			


			Vector3f8 delta_ai_avx;
			delta_ai_avx.setZero();
			Vector3f8 delta_angAcceli_avx;
			delta_angAcceli_avx.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase_avx(
				const Scalarf8 Vj_avx = convert_zero(model->getVolume(0), count);
				compute_Vj_gradW_samephase();

				const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getVelocity(0), count);
				const Vector3f8 omegaj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &m_omega[0], count);

				// Viscosity
				const Scalarf8 density_j_avx = convert_one(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getDensity(0), count);
				const Vector3f8 omegaij = omegai_avx - omegaj_avx;

				// XSPH for angular velocity field
				const Scalarf8 mj_avx = convert_zero(model->getMass(0), count);
				delta_angAcceli_avx -= omegaij * factor_avx * (Vj_avx / density_j_avx) * CubicKernel_AVX::W(xi_avx - xj_avx);

				// difference curl 
				delta_ai_avx += (omegaij % V_gradW) * (nut_density_i * density0_avx);
				delta_angAcceli_avx += ((vi_avx - vj_avx) % V_gradW) * (nut_density_i_intertiaInverse * density0_avx);
			);

			ai[0] += delta_ai_avx.x().reduce();
			ai[1] += delta_ai_avx.y().reduce();
			ai[2] += delta_ai_avx.z().reduce();
			angAcceli[0] += delta_angAcceli_avx.x().reduce();
			angAcceli[1] += delta_angAcceli_avx.y().reduce();
			angAcceli[2] += delta_angAcceli_avx.z().reduce();

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

#else

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

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &vj = m_model->getVelocity(neighborIndex);
				const Vector3r &omegaj = m_omega[neighborIndex];

				// Viscosity
				const Real density_j = m_model->getDensity(neighborIndex);
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

			angAcceli -= 2.0 * m_inertiaInverse * nu_t * omegai;
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < static_cast<int>(numParticles); i++)
		{
			m_omega[i] += dt*m_angularAcceleration[i];
		}
	}
}

#endif

void MicropolarModel_Bender2017::reset()
{
	for (unsigned int i = 0; i < m_model->numParticles(); i++)
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

