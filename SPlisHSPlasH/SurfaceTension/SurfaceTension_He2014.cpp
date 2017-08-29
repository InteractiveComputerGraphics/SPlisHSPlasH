#include "SurfaceTension_He2014.h"
#include <iostream>

using namespace SPH;

SurfaceTension_He2014::SurfaceTension_He2014(FluidModel *model) :
	SurfaceTensionBase(model)
{
	m_color.resize(model->numParticles(), 0.0);
	m_gradC2.resize(model->numParticles(), 0.0);
}

SurfaceTension_He2014::~SurfaceTension_He2014(void)
{
	m_color.clear();
	m_gradC2.clear();
}


void SurfaceTension_He2014::step()
{
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = getSurfaceTension();
	const Real density0 = m_model->getDensity0();

	// Compute color field
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			Real &ci = getColor(i);
			ci = m_model->getMass(i) / m_model->getDensity(i) * m_model->W_zero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Real density_j = m_model->getDensity(neighborIndex);
				ci += m_model->getMass(neighborIndex) / density_j * m_model->W(xi - xj);
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
					ci += m_model->getBoundaryPsi(pid, neighborIndex) / density0 * m_model->W(xi - xj);
				}
			}
		}
	}

	// Compute gradient of color field
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			Vector3r gradC_i;
			gradC_i.setZero();
			const Real &density_i = m_model->getDensity(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Real &density_j = m_model->getDensity(neighborIndex);
				gradC_i += m_model->getMass(neighborIndex) / density_j * getColor(neighborIndex) * m_model->gradW(xi - xj);
			}
			gradC_i *= (1.0 / getColor(i));
			Real &gradC2_i = getGradC2(i);
			gradC2_i = gradC_i.squaredNorm();
		}
	}

	// Compute surface tension force
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			const Real &gradC2_i = getGradC2(i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real &density_i = m_model->getDensity(i);
			const Real factor = 0.25*k / density_i;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Real &gradC2_j = getGradC2(neighborIndex);
				const Real &density_j = m_model->getDensity(neighborIndex);
				ai += factor*m_model->getMass(neighborIndex) / density_j * (gradC2_i + gradC2_j) * m_model->gradW(xi - xj);
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
					ai += factor*m_model->getBoundaryPsi(pid, neighborIndex) / density0 * gradC2_i * m_model->gradW(xi - xj);
				}
			}
		}
	}
}


void SurfaceTension_He2014::reset()
{
}

void SPH::SurfaceTension_He2014::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	auto const& d = m_model->getNeighborhoodSearch()->point_set(0);
	d.sort_field(&m_color[0]);
	d.sort_field(&m_gradC2[0]);
}
