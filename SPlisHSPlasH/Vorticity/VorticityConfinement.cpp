#include "VorticityConfinement.h"
#include <iostream>
#include "../TimeManager.h"
#include "../Simulation.h"

using namespace SPH;
using namespace GenParam;

std::string VorticityConfinement::METHOD_NAME = "Vorticity confinement";
int VorticityConfinement::VORTICITY_COEFFICIENT = -1;

VorticityConfinement::VorticityConfinement(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_vorticityCoeff = static_cast<Real>(0.01);

	m_omega.resize(model->numParticles(), Vector3r::Zero());
	m_normOmega.resize(model->numParticles(), 0.0);

	model->addField({ "angular velocity", METHOD_NAME, FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_omega[i][0]; } });
}

VorticityConfinement::~VorticityConfinement(void)
{
	m_model->removeFieldByName("angular velocity");

	m_omega.clear();
	m_normOmega.clear();
}

void VorticityConfinement::initParameters()
{
	NonPressureForceBase::initParameters();

	VORTICITY_COEFFICIENT = createNumericParameter("vorticity", "Vorticity coefficient", &m_vorticityCoeff);
	setGroup(VORTICITY_COEFFICIENT, "Fluid Model|Vorticity");
	setDescription(VORTICITY_COEFFICIENT, "Coefficient for the vorticity force computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VORTICITY_COEFFICIENT));
	rparam->setMinValue(0.0);
}

void VorticityConfinement::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel *model = m_model;

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			Vector3r &omegai = m_omega[i];
			omegai.setZero();
			const Real density_i = m_model->getDensity(i);
			const Real density_i2 = density_i *density_i;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &vj = m_model->getVelocity(neighborIndex);
				const Real density_j = m_model->getDensity(neighborIndex);
				const Real density_j2 = density_j *density_j;
				const Vector3r gradW = sim->gradW(xi - xj);

				omegai -= m_model->getMass(neighborIndex) / density_i * (vi - vj).cross(gradW);
			)
			Real &normOmegai = m_normOmega[i];
			normOmegai = omegai.norm();
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real density_i = m_model->getDensity(i);
			const Real density_i2 = density_i *density_i;
			const Vector3r &omegai = m_omega[i];

			Vector3r etai;
			etai.setZero();
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = m_model->getDensity(neighborIndex);
				const Vector3r gradW = sim->gradW(xi - xj);
				Real &normOmegaj = m_normOmega[neighborIndex];
				etai += m_model->getMass(neighborIndex) / density_i * normOmegaj * gradW;
			)

			etai.normalize();
			ai += m_vorticityCoeff * etai.cross(omegai);
		}
	}
}


void VorticityConfinement::reset()
{
}

void VorticityConfinement::performNeighborhoodSearchSort()
{
}

