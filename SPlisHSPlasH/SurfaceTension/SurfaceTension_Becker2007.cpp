#include "SurfaceTension_Becker2007.h"
#include <iostream>

using namespace SPH;

SurfaceTension_Becker2007::SurfaceTension_Becker2007(FluidModel *model) :
	SurfaceTensionBase(model)
{
}

SurfaceTension_Becker2007::~SurfaceTension_Becker2007(void)
{
}

void SurfaceTension_Becker2007::step()
{
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = getSurfaceTension();
	const Real diameter = 2.0 * m_model->getParticleRadius();
	const Real diameter2 = diameter*diameter;

	// Compute forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			Vector3r &ai = m_model->getAcceleration(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Vector3r xixj = xi - xj;
				const Real r2 = xixj.dot(xixj);
				if (r2 > diameter2)
					ai -= k / m_model->getMass(i) * m_model->getMass(neighborIndex) * (xi - xj) * m_model->W(xi - xj);
				else
					ai -= k / m_model->getMass(i) * m_model->getMass(neighborIndex) * (xi - xj) * m_model->W(Vector3r(diameter, 0.0, 0.0));
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
			{
				for (unsigned int j = 0; j < m_model->numberOfNeighbors(pid, i); j++)
				{
					const unsigned int neighborIndex = m_model->getNeighbor(pid, i, j);
					const Vector3r &xj = m_model->getPosition(pid, neighborIndex);
					const Vector3r xixj = xi - xj;
					const Real r2 = xixj.dot(xixj);
					if (r2 > diameter2)
						ai -= k / m_model->getMass(i) * m_model->getBoundaryPsi(pid, neighborIndex) * (xi - xj) * m_model->W(xi - xj);
					else
						ai -= k / m_model->getMass(i) * m_model->getBoundaryPsi(pid, neighborIndex) * (xi - xj) * m_model->W(Vector3r(diameter, 0.0, 0.0));
				}
			}
		}
	}
}


void SurfaceTension_Becker2007::reset()
{
}

