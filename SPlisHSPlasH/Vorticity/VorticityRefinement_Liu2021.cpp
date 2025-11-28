#include "VorticityRefinement_Liu2021.h"
#include <iostream>
#include "../TimeManager.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace GenParam;

std::string VorticityRefinement_Liu2021::METHOD_NAME = "Liu et al. 2021";
int VorticityRefinement_Liu2021::VORTICITY_COEFFICIENT = -1;
int VorticityRefinement_Liu2021::KINEMATIC_VISCOSITY = -1;


VorticityRefinement_Liu2021::VorticityRefinement_Liu2021(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_vorticityCoeff = static_cast<Real>(1.0);
	m_viscosityCoefficient = static_cast<Real>(0.05);

	m_vorticity.resize(model->numParticles(), Vector3r::Zero());
	m_dissipatedVorticity.resize(model->numParticles(), Vector3r::Zero());
	m_streamFunction.resize(model->numParticles(), Vector3r::Zero());

	model->addField({ "vorticity", METHOD_NAME, FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity[i][0]; }, true });
	model->addField({ "dissipated vorticity", METHOD_NAME, FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_dissipatedVorticity[i][0]; }, true });
	model->addField({ "stream function", METHOD_NAME, FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_streamFunction[i][0]; }, true });
}

VorticityRefinement_Liu2021::~VorticityRefinement_Liu2021(void)
{
	m_model->removeFieldByName("vorticity");
	m_model->removeFieldByName("dissipated vorticity");
	m_model->removeFieldByName("stream function");

	m_vorticity.clear();
	m_dissipatedVorticity.clear();
	m_streamFunction.clear();
}

void VorticityRefinement_Liu2021::initParameters()
{
	NonPressureForceBase::initParameters();

	VORTICITY_COEFFICIENT = createNumericParameter("vorticity", "Vorticity coefficient", &m_vorticityCoeff);
	setGroup(VORTICITY_COEFFICIENT, "Fluid Model|Vorticity");
	setDescription(VORTICITY_COEFFICIENT, "Coefficient for the vorticity force computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VORTICITY_COEFFICIENT));
	rparam->setMinValue(0.0);

	KINEMATIC_VISCOSITY = createNumericParameter("kinetmaticViscosity", "Viscosity coefficient", &m_viscosityCoefficient);
	setGroup(KINEMATIC_VISCOSITY, "Fluid Model|Vorticity");
	setDescription(KINEMATIC_VISCOSITY, "Coefficient for the viscosity force computation");
	rparam = static_cast<RealParameter*>(getParameter(KINEMATIC_VISCOSITY));
	rparam->setMinValue(0.0);
}

void VorticityRefinement_Liu2021::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel *model = m_model;
	const Real density0 = model->getDensity0();

	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real invDt = static_cast<Real>(1.0) / dt;

	const Real alpha = m_vorticityCoeff;

	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;

	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 8.0;

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// compute vorticity through linear field 
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			Vector3r &vorticity_i = m_vorticity[i];
			Matrix3r grad_v = Matrix3r::Zero();
			Vector3r vorticityLaplacian = Vector3r::Zero();
			Vector3r vorticity_linear_field = Vector3r::Zero();
	
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &vj = m_model->getVelocity(neighborIndex);
				const Real density_j = m_model->getDensity(neighborIndex);
				const Vector3r xij = xi - xj;
				const Vector3r vij = vi - vj;
				const Vector3r gradW = sim->gradW(xij);
				const Real Vj = m_model->getMass(neighborIndex) / density_j;
				const Vector3r& vorticity_j = m_vorticity[neighborIndex];

				// compute vorticity through linear field 
				vorticity_linear_field += Vj * vij.cross(gradW);
				// compute grad v
				grad_v -= Vj * vij * gradW.transpose();
				// compute Laplacian of vorticity
				vorticityLaplacian += d * Vj * (vorticity_i - vorticity_j).dot(xij) / (xij.squaredNorm() + 0.01 * h2) * gradW;
			);

			const Vector3r strechingTerm = grad_v * vorticity_i;
			const Vector3r vorticityPrediction = vorticity_i + dt * (strechingTerm + m_viscosityCoefficient * vorticityLaplacian);

			Vector3r& dissipatedVorticity_i = m_dissipatedVorticity[i];
			dissipatedVorticity_i = vorticityPrediction - vorticity_linear_field;
		}

		//////////////////////////////////////////////////////////////////////////
		// compute stream function
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < static_cast<int>(numParticles); i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			Vector3r& streamFunction_i = m_streamFunction[i];
			streamFunction_i.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			const Real epsilon = static_cast<Real>(0.1) * h; 
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = m_model->getDensity(neighborIndex);
				const Vector3r& dissipatedVorticity_j = m_dissipatedVorticity[neighborIndex];
				const Real Vj = m_model->getMass(neighborIndex) / density_j;
				streamFunction_i += dissipatedVorticity_j * Vj / ((xi - xj).norm() + epsilon);		// avoid a division by zero using epsilon in the denominator
			);

			streamFunction_i *= static_cast<Real>(0.25 / M_PI);
		}

		//////////////////////////////////////////////////////////////////////////
		// compute refinement of linear velocity
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < static_cast<int>(numParticles); i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			const Vector3r& streamFunction_i = m_streamFunction[i];
			Vector3r& vi = m_model->getVelocity(i);

			Vector3r delta_vi = Vector3r::Zero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = m_model->getDensity(neighborIndex);
				const Vector3r& streamFunction_j = m_streamFunction[neighborIndex];
				const Vector3r gradW = sim->gradW(xi - xj);
				const Real Vj = m_model->getMass(neighborIndex) / density_j;

				delta_vi += Vj * (streamFunction_i - streamFunction_j).cross(gradW);
			);

			if (m_model->getParticleState(i) == ParticleState::Active)
				vi += m_vorticityCoeff * delta_vi;
		}

		//////////////////////////////////////////////////////////////////////////
		// store the vorticity for the next time step
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static)  
		for (int i = 0; i < static_cast<int>(numParticles); i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			const Vector3r& vi = m_model->getVelocity(i);
			Vector3r& vorticity_i = m_vorticity[i];
			vorticity_i.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////

			forall_fluid_neighbors_in_same_phase(
				const Vector3r & vj = m_model->getVelocity(neighborIndex);
				const Real density_j = m_model->getDensity(neighborIndex);
				const Vector3r gradW = sim->gradW(xi - xj);

				vorticity_i += m_model->getMass(neighborIndex) / density_j * (vi - vj).cross(gradW);
			);
		}
	}
}

void VorticityRefinement_Liu2021::reset()
{
	for (unsigned int i = 0; i < m_model->numParticles(); i++)
	{
		m_vorticity[i].setZero();
		m_dissipatedVorticity[i].setZero();
		m_streamFunction[i].setZero();
	}
}

void SPH::VorticityRefinement_Liu2021::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_vorticity[0]);
}

