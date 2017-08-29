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
	const unsigned int numParticles = m_model->numActiveParticles();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = (1.0 / h);

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

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Vector3r &vj = m_model->getVelocity(0, neighborIndex);

				// Viscosity
				const Real density_j = m_model->getDensity(neighborIndex);
				ai -= invH * m_viscosity * (m_model->getMass(neighborIndex) / density_j) * (vi - vj) * m_model->W(xi - xj);
			}

			////////////////////////////////////////////////////////////////////////////
			//// Boundary
			////////////////////////////////////////////////////////////////////////////
			//for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
			//{
			//	for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
			//	{
			//		const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
			//		const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
			//		const Vector3r &vj = m_model->getVelocity(pid, neighborIndex);
			//		ai -= invH * viscosity * (m_model->getBoundaryPsi(pid, neighborIndex) / density_i) * (vi)* m_model->W(xi - xj);
			//	}
			//}
		}
	}
}


void Viscosity_XSPH::reset()
{
}

