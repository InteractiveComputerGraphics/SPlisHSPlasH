#include "Viscosity_XSPH.h"
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;

Viscosity_XSPH::Viscosity_XSPH(FluidModel *model) :
	ViscosityBase(model)
{
}

Viscosity_XSPH::~Viscosity_XSPH(void)
{
}

void Viscosity_XSPH::step()
{
	const unsigned int numParticles = m_model->numParticles();

	const Real viscosity = m_model->getViscosity();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Compute viscosity forces (XSPH)
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			const Vector3r &vi = m_model->getVelocity(0, i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real density_i = m_model->getDensity(i);
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
				const unsigned int &neighborIndex = particleId.point_id;
				const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);
				const Vector3r &vj = m_model->getVelocity(particleId.point_set_id, neighborIndex);

				if (particleId.point_set_id == 0)		// Test if fluid particle
				{
					// Viscosity
					const Real density_j = m_model->getDensity(neighborIndex);
					ai -= (1.0/h) * viscosity * (m_model->getMass(neighborIndex) / density_j) * (vi - vj) * m_model->W(xi - xj);

				}
// 				else 
// 				{
// 					ai -= (1.0/h) * viscosity * (m_model->getBoundaryPsi(pid, neighborIndex) / density_i) * (vi)* m_model->W(xi - xj);
// 				}
			}
		}
	}
}


void Viscosity_XSPH::reset()
{
}

