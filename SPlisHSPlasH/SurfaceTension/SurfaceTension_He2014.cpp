#include "SurfaceTension_He2014.h"
#include <iostream>
#include "../Simulation.h"

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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const Real density0 = m_model->getValue<Real>(FluidModel::DENSITY0);
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel *model = m_model;

	// Compute color field
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Real &ci = getColor(i);
			ci = m_model->getMass(i) / m_model->getDensity(i) * sim->W_zero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = m_model->getDensity(neighborIndex);
				ci += m_model->getMass(neighborIndex) / density_j * sim->W(xi - xj);
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			forall_boundary_neighbors(
					ci += bm_neighbor->getBoundaryPsi(neighborIndex) / density0 * sim->W(xi - xj);
			)
		}
	}

	// Compute gradient of color field
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r gradC_i;
			gradC_i.setZero();
			const Real &density_i = m_model->getDensity(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real &density_j = m_model->getDensity(neighborIndex);
				gradC_i += m_model->getMass(neighborIndex) / density_j * getColor(neighborIndex) * sim->gradW(xi - xj);
			)
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
			const Vector3r &xi = m_model->getPosition(i);
			const Real &gradC2_i = getGradC2(i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real &density_i = m_model->getDensity(i);
			const Real factor = 0.25*k / density_i;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real &gradC2_j = getGradC2(neighborIndex);
				const Real &density_j = m_model->getDensity(neighborIndex);
				ai += factor*m_model->getMass(neighborIndex) / density_j * (gradC2_i + gradC2_j) * sim->gradW(xi - xj);
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			forall_boundary_neighbors(
					ai += factor*bm_neighbor->getBoundaryPsi(neighborIndex) / density0 * gradC2_i * sim->gradW(xi - xj);
			)
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

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_color[0]);
	d.sort_field(&m_gradC2[0]);
}
