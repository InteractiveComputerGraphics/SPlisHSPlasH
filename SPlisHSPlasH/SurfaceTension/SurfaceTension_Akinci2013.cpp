#include "SurfaceTension_Akinci2013.h"
#include <iostream>
#include "../Simulation.h"

using namespace SPH;

SurfaceTension_Akinci2013::SurfaceTension_Akinci2013(FluidModel *model) :
	SurfaceTensionBase(model)
{
	m_normals.resize(model->numParticles(), Vector3r::Zero());

	model->addField({ "normal", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_normals[i][0]; } });
}

SurfaceTension_Akinci2013::~SurfaceTension_Akinci2013(void)
{
	m_model->removeFieldByName("normal");
	m_normals.clear();
}


void SurfaceTension_Akinci2013::computeNormals()
{
	Simulation *sim = Simulation::getCurrent();
	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel *model = m_model;

	// Compute normals
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r &ni = getNormal(i);
			ni.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = m_model->getDensity(neighborIndex);
				ni += m_model->getMass(neighborIndex) / density_j * sim->gradW(xi - xj);
			)
			ni = supportRadius*ni;
		}
	}

}

void SurfaceTension_Akinci2013::step()
{
	Simulation *sim = Simulation::getCurrent();
	const Real density0 = m_model->getDensity0();
	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel *model = m_model;

	computeNormals();

	// Compute forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &ni = getNormal(i);
			const Real &rhoi = m_model->getDensity(i);
			Vector3r &ai = m_model->getAcceleration(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real &rhoj = m_model->getDensity(neighborIndex);
				const Real K_ij = static_cast<Real>(2.0)*density0 / (rhoi + rhoj);

				Vector3r accel;
				accel.setZero();

				// Cohesion force
				Vector3r xixj = (xi - xj);
				const Real length2 = xixj.squaredNorm();
				if (length2 > 1.0e-9)
				{
					xixj = (static_cast<Real>(1.0) / sqrt(length2)) * xixj;
					accel -= k * m_model->getMass(neighborIndex) * xixj * CohesionKernel::W(xi - xj);
				}

				// Curvature
				const Vector3r &nj = getNormal(neighborIndex);
				accel -= k * supportRadius* (ni - nj);

				ai += K_ij * accel;
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			forall_boundary_neighbors(
					// adhesion force					
					Vector3r xixj = (xi - xj);
					const Real length2 = xixj.squaredNorm();
					if (length2 > 1.0e-9)
					{
						xixj = ((Real) 1.0 / sqrt(length2)) * xixj;
						ai -= k * density0 * bm_neighbor->getVolume(neighborIndex) * xixj * AdhesionKernel::W(xi - xj);
					}
			)
		}
	}
}


void SurfaceTension_Akinci2013::reset()
{
}

void SurfaceTension_Akinci2013::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_normals[0]);
}

