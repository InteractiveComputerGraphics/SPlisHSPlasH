#include "Viscosity_Standard.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../Simulation.h"

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
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 6.0;

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
				const Vector3r xixj = xi - xj;
				ai += d * m_viscosity * (fm_neighbor->getMass(neighborIndex) / density_j) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * sim->gradW(xi - xj);
			)

			////////////////////////////////////////////////////////////////////////////
			//// Boundary
			////////////////////////////////////////////////////////////////////////////
			//for (unsigned int pid = nFluids; pid < m_model->numberOfPointSets(); pid++)
			//{
			//  BoundaryModel *bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
			//	for (unsigned int j = 0; j < m_model->numberOfNeighbors(fluidModelIndex, pid, i); j++)
			//	{
			//		const unsigned int neighborIndex = m_model->getNeighbor(fluidModelIndex, pid, i, j);
			//		const Vector3r &xj = bm_neighbor->getPosition(neighborIndex);
			//		const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
			//		const Vector3r xixj = xi - xj;
			//		ai += d * viscosity * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi) * (xixj.dot(sim->gradW(xi - xj))) / (xixj.squaredNorm() + 0.01*h2);
			//	}
			//}
		}
	}
}


void Viscosity_Standard::reset()
{
}

