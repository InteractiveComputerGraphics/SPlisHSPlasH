#include "XSPH.h"
#include "TimeManager.h"
#include "Simulation.h"
#include "BoundaryModel_Akinci2012.h"
#include "BoundaryModel_Koschier2017.h"
#include "BoundaryModel_Bender2019.h"

using namespace SPH;
using namespace GenParam;

int XSPH::FLUID_COEFFICIENT = -1;
int XSPH::BOUNDARY_COEFFICIENT = -1;


XSPH::XSPH(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_fluidCoefficient = 0.0;
	m_boundaryCoefficient = 0.0;
}

XSPH::~XSPH(void)
{
}

void XSPH::initParameters()
{
	NonPressureForceBase::initParameters();

	FLUID_COEFFICIENT = createNumericParameter("xsph", "Fluid coefficient", &m_fluidCoefficient);
	setGroup(FLUID_COEFFICIENT, "Fluid Model|XSPH");
	setDescription(FLUID_COEFFICIENT, "Coefficient for the XSPH velocity filter in the fluid");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(FLUID_COEFFICIENT));
	rparam->setMinValue(0.0);

	BOUNDARY_COEFFICIENT = createNumericParameter("xsphBoundary", "Boundary coefficient", &m_boundaryCoefficient);
	setGroup(BOUNDARY_COEFFICIENT, "Fluid Model|XSPH");
	setDescription(BOUNDARY_COEFFICIENT, "Coefficient for the XSPH velocity filter at the boundary");
	rparam = static_cast<RealParameter*>(getParameter(BOUNDARY_COEFFICIENT));
	rparam->setMinValue(0.0);
}

#ifdef USE_AVX

void XSPH::step()
{
	if ((m_fluidCoefficient == 0.0) && (m_boundaryCoefficient == 0.0))
		return;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real density0 = m_model->getValue<Real>(FluidModel::DENSITY0);

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

			const Vector3f8 xi_avx(xi);
			const Vector3f8 vi_avx(vi);
			const Scalarf8 mi_avx(m_model->getMass(i));
			const Scalarf8 density_i_avx(density_i);
			Vector3f8 delta_ai;
			delta_ai.setZero();
			const Scalarf8 dvisc(invH * m_fluidCoefficient);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			////////////////////////////////////////////////////////////////////////
			if (m_fluidCoefficient != 0.0)
			{
				forall_fluid_neighbors_avx(
					const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getVelocity(0), count);
					const Scalarf8 mj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getMass(0), count);

					// Viscosity
					const Scalarf8 density_j_avx = convert_one(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getDensity(0), count);
					const Vector3f8 xixj = xi_avx - xj_avx;
					delta_ai -= (vi_avx - vj_avx) * dvisc * (mj_avx / density_j_avx) * CubicKernel_AVX::W(xixj);
				);
				ai += delta_ai.reduce();
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (m_boundaryCoefficient != 0.0)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
						const Vector3r a = -invH * m_boundaryCoefficient * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi-vj)* sim->W(xi - xj);
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						Vector3r vj;
						bm_neighbor->getPointVelocity(xi, vj);
						const Vector3r a = -invH * m_boundaryCoefficient * (density0 / density_i) * (vi-vj)* rho;
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						Vector3r vj;
						bm_neighbor->getPointVelocity(xj, vj);
						const Vector3r a = -invH * m_boundaryCoefficient * (density0 * Vj / density_i) * (vi-vj)* sim->W(xi - xj);
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
			}
		}
	}
}

#else

void XSPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real density0 = m_model->getValue<Real>(FluidModel::DENSITY0);

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
				ai -= invH * m_fluidCoefficient * (fm_neighbor->getMass(neighborIndex) / density_j) * (vi - vj) * sim->W(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (m_boundaryCoefficient != 0.0)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
						const Vector3r a = -invH * m_boundaryCoefficient * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi-vj)* sim->W(xi - xj);
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						Vector3r vj;
						bm_neighbor->getPointVelocity(xi, vj);
						const Vector3r a = -invH * m_boundaryCoefficient * (density0 / density_i) * (vi-vj)* rho;
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						Vector3r vj;
						bm_neighbor->getPointVelocity(xj, vj);
						const Vector3r a = -invH * m_boundaryCoefficient * (density0 * Vj / density_i) * (vi-vj)* sim->W(xi - xj);
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
			}
		}
	}
}

#endif

void XSPH::reset()
{
}

