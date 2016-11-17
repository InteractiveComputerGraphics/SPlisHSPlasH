#include "Viscosity_Standard.h"
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;

Viscosity_Standard::Viscosity_Standard(FluidModel *model) :
	ViscosityBase(model)
{
}

Viscosity_Standard::~Viscosity_Standard(void)
{
}

void Viscosity_Standard::step()
{
	const unsigned int numParticles = m_model->numParticles();
	const Real h = m_model->getSupportRadius();
	const Real h2 = h*h;
	const Real viscosity = m_model->getViscosity();

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
					const Vector3r xixj = xi - xj;
					ai += 2.0 * viscosity * (m_model->getMass(neighborIndex) / density_j) * (vi - vj) * (xixj.dot(m_model->gradW(xi - xj)))/(xixj.squaredNorm() + 0.01*h2);
				}
// 				else 
// 				{
//					ai += 2.0 * viscosity * (model.getBoundaryPsi(particleId.point_set_id, neighborIndex) / density_i) * (vi) * (xixj.dot(m_model->gradW(xi - xj))) / (xixj.squaredNorm() + 0.01*h2);
// 				}
			}
		}
	}
}


void Viscosity_Standard::reset()
{
}

