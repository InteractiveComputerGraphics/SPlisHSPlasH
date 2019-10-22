#include "Viscosity_XSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;
using namespace GenParam;

int Viscosity_XSPH::VISCOSITY_COEFFICIENT_BOUNDARY = -1;


Viscosity_XSPH::Viscosity_XSPH(FluidModel *model) :
	ViscosityBase(model)
{
	m_boundaryViscosity = 0.0;
}

Viscosity_XSPH::~Viscosity_XSPH(void)
{
}

void Viscosity_XSPH::initParameters()
{
	ViscosityBase::initParameters();

	VISCOSITY_COEFFICIENT_BOUNDARY = createNumericParameter("viscosityBoundary", "Viscosity coefficient (Boundary)", &m_boundaryViscosity);
	setGroup(VISCOSITY_COEFFICIENT_BOUNDARY, "Viscosity");
	setDescription(VISCOSITY_COEFFICIENT_BOUNDARY, "Coefficient for the viscosity force computation at the boundary.");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT_BOUNDARY));
	rparam->setMinValue(0.0);
}

void Viscosity_XSPH::step()
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
				ai -= invH * m_viscosity * (fm_neighbor->getMass(neighborIndex) / density_j) * (vi - vj) * sim->W(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (m_boundaryViscosity != 0.0)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
						const Vector3r a = -invH * m_boundaryViscosity * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi-vj)* sim->W(xi - xj);
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						Vector3r vj;
						bm_neighbor->getPointVelocity(xi, vj);
						const Vector3r a = -invH * m_boundaryViscosity * (density0 / density_i) * (vi-vj)* rho;
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						Vector3r vj;
						bm_neighbor->getPointVelocity(xj, vj);
						const Vector3r a = -invH * m_boundaryViscosity * (density0 * Vj / density_i) * (vi-vj)* sim->W(xi - xj);
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
			}
		}
	}
}


void Viscosity_XSPH::reset()
{
}

