#include "SurfaceTension_Becker2007.h"
#include <iostream>
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;
using namespace GenParam;

std::string SurfaceTension_Becker2007::METHOD_NAME = "Becker & Teschner 2007";
int SurfaceTension_Becker2007::SURFACE_TENSION = -1;
int SurfaceTension_Becker2007::SURFACE_TENSION_BOUNDARY = -1;


SurfaceTension_Becker2007::SurfaceTension_Becker2007(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_surfaceTension = static_cast<Real>(0.05);
	m_surfaceTensionBoundary = static_cast<Real>(0.01);
}

SurfaceTension_Becker2007::~SurfaceTension_Becker2007(void)
{
}

void SurfaceTension_Becker2007::initParameters()
{
	NonPressureForceBase::initParameters();

	SURFACE_TENSION = createNumericParameter("surfaceTension", "Surface tension coefficient", &m_surfaceTension);
	setGroup(SURFACE_TENSION, "Fluid Model|Surface tension");
	setDescription(SURFACE_TENSION, "Coefficient for the surface tension computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(SURFACE_TENSION));
	rparam->setMinValue(0.0);

	SURFACE_TENSION_BOUNDARY = createNumericParameter("surfaceTensionBoundary", "Boundary surface tension coefficient", &m_surfaceTensionBoundary);
	setGroup(SURFACE_TENSION_BOUNDARY, "Fluid Model|Surface tension");
	setDescription(SURFACE_TENSION_BOUNDARY, "Coefficient for the surface tension computation at the boundary");
	rparam = static_cast<RealParameter*>(getParameter(SURFACE_TENSION_BOUNDARY));
	rparam->setMinValue(0.0);
}


void SurfaceTension_Becker2007::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const Real kb = m_surfaceTensionBoundary;
	const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
	const Real diameter = static_cast<Real>(2.0) * radius;
	const Real diameter2 = diameter*diameter;
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = m_model->getDensity0();
	FluidModel *model = m_model;
	const Real h = sim->getSupportRadius();

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
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r xixj = xi - xj;
					const Real r2 = xixj.dot(xixj);
						Vector3r a;
					if (r2 > diameter2)
							a = -kb / m_model->getMass(i) * density0 * bm_neighbor->getVolume(neighborIndex) * (xi - xj) * sim->W(xi - xj);
					else
							a = -kb / m_model->getMass(i) * density0 * bm_neighbor->getVolume(neighborIndex) * (xi - xj) * sim->W(Vector3r(diameter, 0.0, 0.0));
						ai += a;
						bm_neighbor->addForce(xj, -model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					const Vector3r xixj = xi - xj;
					const Real r2 = xixj.dot(xixj);
					Vector3r a;
					if (r2 > diameter2)
						a = -kb / m_model->getMass(i) * rho * density0 * (xi - xj) * rho;
 					else
 						a = -kb / m_model->getMass(i) * rho * density0 * (xi - xj) * sim->W(Vector3r(diameter, 0.0, 0.0)) / sim->W(xi - xj);
					ai += a;
					bm_neighbor->addForce(xj, -model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r xixj = xi - xj;
					const Real r2 = xixj.dot(xixj);
					Vector3r a;
					a.setZero();
					if (r2 > diameter2)
						a -= kb / m_model->getMass(i) * Vj * density0 * (xi - xj) * sim->W(xi - xj);
					else
						a -= kb / m_model->getMass(i) * Vj * density0 * (xi - xj) * sim->W(Vector3r(diameter, 0.0, 0.0));
					ai += a;
					bm_neighbor->addForce(xj, -model->getMass(i) * a);
				);
			}
		}
	}
}


void SurfaceTension_Becker2007::reset()
{
}

