#include "SurfaceTension_Akinci2013.h"
#include <iostream>

using namespace SPH;

SurfaceTension_Akinci2013::SurfaceTension_Akinci2013(FluidModel *model) :
	SurfaceTensionBase(model)
{
	m_normals.resize(model->numParticles(), SPH::Vector3r::Zero());
}

SurfaceTension_Akinci2013::~SurfaceTension_Akinci2013(void)
{
	m_normals.clear();
}


void SurfaceTension_Akinci2013::computeNormals()
{
	const Real supportRadius = m_model->getSupportRadius();
	const unsigned int numParticles = m_model->numParticles();

	// Compute normals
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			Vector3r &ni = getNormal(i);
			ni.setZero();
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
				const unsigned int &neighborIndex = particleId.point_id;

				if (particleId.point_set_id == 0)
				{
					const Vector3r &xj = m_model->getPosition(0, neighborIndex);	
					const Real density_j = m_model->getDensity(neighborIndex);
					ni += m_model->getMass(neighborIndex) / density_j * m_model->gradW(xi - xj);
				}
			}
			ni = supportRadius*ni;
		}
	}

}

void SurfaceTension_Akinci2013::step()
{
	const Real density0 = m_model->getDensity0();
	const Real supportRadius = m_model->getSupportRadius();
	const unsigned int numParticles = m_model->numParticles();
	const Real k = m_model->getSurfaceTension();

	computeNormals();

	// Compute forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			const Vector3r &ni = getNormal(i);
			const Real &rhoi = m_model->getDensity(i);
			Vector3r &ai = m_model->getAcceleration(i);
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
				const unsigned int &neighborIndex = particleId.point_id;
				const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

				if (particleId.point_set_id == 0)
				{
					const Real &rhoj = m_model->getDensity(neighborIndex);
					const Real K_ij = 2.0*density0 / (rhoi + rhoj);

					Vector3r accel;
					accel.setZero();

					// Cohesion force
					Vector3r xixj = (xi - xj);
					const Real length2 = xixj.squaredNorm();
					if (length2 > 1.0e-9)
					{
						xixj = ((Real) 1.0 / sqrt(length2)) * xixj;
						accel -= k * m_model->getMass(neighborIndex) * xixj * CohesionKernel::W(xi - xj);
					}

					// Curvature
					const Vector3r &nj = getNormal(neighborIndex);
					accel -= k * supportRadius* (ni - nj);

					ai += K_ij * accel;
				}
				else
				{
					// adhesion force					
					Vector3r xixj = (xi - xj);
					const Real length2 = xixj.squaredNorm();
					if (length2 > 1.0e-9)
					{
						xixj = ((Real) 1.0 / sqrt(length2)) * xixj;
						ai -= k * m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) * xixj * AdhesionKernel::W(xi - xj);
					}
				}
			}
		}
	}
}


void SurfaceTension_Akinci2013::reset()
{
}

void SurfaceTension_Akinci2013::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_normals[0]);
}

