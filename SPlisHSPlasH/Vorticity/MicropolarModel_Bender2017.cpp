#include "MicropolarModel_Bender2017.h"
#include <iostream>
#include "../TimeManager.h"

using namespace SPH;

MicropolarModel_Bender2017::MicropolarModel_Bender2017(FluidModel *model) :
	VorticityBase(model)
{
	m_omega.resize(model->numParticles(), Vector3r::Zero());
	m_angularAcceleration.resize(model->numParticles(), Vector3r::Zero());
	m_inertiaInverse = 0.5;
	m_viscosityOmega = 0.1;
}

MicropolarModel_Bender2017::~MicropolarModel_Bender2017(void)
{
	m_omega.clear();
	m_angularAcceleration.clear();
}


void MicropolarModel_Bender2017::step()
{
	const unsigned int numParticles = m_model->numActiveParticles();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = 1.0 / h;

	const Real nu_t = m_vorticityCoeff;
	const Real zeta = m_viscosityOmega;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			const Vector3r &vi = m_model->getVelocity(0, i);
			const Vector3r &omegai = m_omega[i];
			Vector3r &ai = m_model->getAcceleration(i);
			Vector3r &angAcceli = m_angularAcceleration[i];
			angAcceli.setZero();
			const Real density_i = m_model->getDensity(i);
			const Real density_i2 = density_i *density_i;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Vector3r &vj = m_model->getVelocity(0, neighborIndex);
				const Vector3r &omegaj = m_omega[neighborIndex];

				// Viscosity
				const Real density_j = m_model->getDensity(neighborIndex);
				const Real density_j2 = density_j *density_j;
				const Vector3r xij = xi - xj;
				const Vector3r omegaij = omegai - omegaj;
				const Vector3r gradW = m_model->gradW(xij);

				// XSPH for angular velocity field
				angAcceli -= invH * m_inertiaInverse * zeta * (m_model->getMass(neighborIndex) / density_j) * omegaij * m_model->W(xij);

				// symmetric curl 
				ai -= nu_t * density_i * m_model->getMass(neighborIndex) * ((omegai / density_i2 + omegaj / density_j2).cross(gradW));
				angAcceli -= nu_t * density_i * m_inertiaInverse * (m_model->getMass(neighborIndex) * (vi / density_i2 + vj / density_j2).cross(gradW));
			}
			angAcceli -= 2.0 * m_inertiaInverse * nu_t * omegai;
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_omega[i] += h*m_angularAcceleration[i];
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

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_omega[0]);
}

