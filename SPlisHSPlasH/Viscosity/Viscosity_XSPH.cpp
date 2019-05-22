#include "Viscosity_XSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../Simulation.h"

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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int numParticles = m_model->numActiveParticles();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = (static_cast<Real>(1.0) / h);

	// Compute viscosity forces (XSPH)
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real density_i = m_model->getDensity(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r &vj = fm_neighbor->getVelocity(neighborIndex);

				// Viscosity
				const Real density_j = fm_neighbor->getDensity(neighborIndex);
				ai -= invH * m_viscosity * (fm_neighbor->getMass(neighborIndex) / density_j) * (vi - vj) * sim->W(xi - xj);
			);

			////////////////////////////////////////////////////////////////////////////
			//// Boundary
			////////////////////////////////////////////////////////////////////////////
			//for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
			//{
			//  BoundaryModel *bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
			//	for (unsigned int j = 0; j < sim->numberOfNeighbors(pid, i); j++)
			//	{
			//		const unsigned int neighborIndex = sim->getNeighbor(pid, i, j);
			//		const Vector3r &xj = bm_neighbor->getPosition(neighborIndex);
			//		const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
			//		ai -= invH * viscosity * bm_neighbor->getVolume(neighborIndex) * (vi)* sim->W(xi - xj);
			//	}
			//}
		}
	}
}


void Viscosity_XSPH::reset()
{
}

