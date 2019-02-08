#include "SurfaceTension_Becker2007.h"
#include <iostream>
#include "../Simulation.h"

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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
	const Real diameter = static_cast<Real>(2.0) * radius;
	const Real diameter2 = diameter*diameter;
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel *model = m_model;
	const Real density0 = model->getDensity0();

	// Compute forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r &ai = m_model->getAcceleration(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r xixj = xi - xj;
				const Real r2 = xixj.dot(xixj);
				if (r2 > diameter2)
					ai -= k / m_model->getMass(i) * m_model->getMass(neighborIndex) * (xi - xj) * sim->W(xi - xj);
				else
					ai -= k / m_model->getMass(i) * m_model->getMass(neighborIndex) * (xi - xj) * sim->W(Vector3r(diameter, 0.0, 0.0));
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			forall_boundary_neighbors(
				const Vector3r xixj = xi - xj;
				const Real r2 = xixj.dot(xixj);
				if (r2 > diameter2)
					ai -= k / m_model->getMass(i) * density0 * bm_neighbor->getVolume(neighborIndex) * (xi - xj) * sim->W(xi - xj);
				else
					ai -= k / m_model->getMass(i) * density0 * bm_neighbor->getVolume(neighborIndex) * (xi - xj) * sim->W(Vector3r(diameter, 0.0, 0.0));
			)
		}
	}
}


void SurfaceTension_Becker2007::reset()
{
}

