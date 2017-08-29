#include "VorticityConfinement.h"
#include <iostream>
#include "../TimeManager.h"

using namespace SPH;

VorticityConfinement::VorticityConfinement(FluidModel *model) :
	VorticityBase(model)
{
	m_omega.resize(model->numParticles(), Vector3r::Zero());
	m_normOmega.resize(model->numParticles(), 0.0);
}

VorticityConfinement::~VorticityConfinement(void)
{
	m_omega.clear();
	m_normOmega.clear();
}


void VorticityConfinement::step()
{
	const unsigned int numParticles = m_model->numActiveParticles();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			const Vector3r &vi = m_model->getVelocity(0, i);
			Vector3r &omegai = m_omega[i];
			omegai.setZero();
			Vector3r &ai = m_model->getAcceleration(i);
			const Real density_i = m_model->getDensity(i);
			const Real density_i2 = density_i *density_i;

			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Vector3r &vj = m_model->getVelocity(0, neighborIndex);
				const Real density_j = m_model->getDensity(neighborIndex);
				const Real density_j2 = density_j *density_j;
				const Vector3r gradW = m_model->gradW(xi - xj);

				omegai -= m_model->getMass(neighborIndex) / density_i * (vi - vj).cross(gradW);
			}
			Real &normOmegai = m_normOmega[i];
			normOmegai = omegai.norm();
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			const Vector3r &vi = m_model->getVelocity(0, i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real density_i = m_model->getDensity(i);
			const Real density_i2 = density_i *density_i;
			const Vector3r &omegai = m_omega[i];

			Vector3r etai;
			etai.setZero();
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Real density_j = m_model->getDensity(neighborIndex);
				const Vector3r gradW = m_model->gradW(xi - xj);
				Real &normOmegaj = m_normOmega[neighborIndex];
				etai += m_model->getMass(neighborIndex) / density_i * normOmegaj * gradW;
			}

			etai.normalize();
			ai += m_vorticityCoeff * etai.cross(omegai);
		}
	}
}


void VorticityConfinement::reset()
{
}

void VorticityConfinement::performNeighborhoodSearchSort()
{
}

