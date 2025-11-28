#include "DragForce_Macklin2014.h"
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;
using namespace GenParam;

std::string DragForce_Macklin2014::METHOD_NAME = "Macklin et al. 2014";
int DragForce_Macklin2014::DRAG_COEFFICIENT = -1;

DragForce_Macklin2014::DragForce_Macklin2014(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_dragCoefficient = static_cast<Real>(0.01);
}

DragForce_Macklin2014::~DragForce_Macklin2014(void)
{
}

void DragForce_Macklin2014::initParameters()
{
	NonPressureForceBase::initParameters();

	DRAG_COEFFICIENT = createNumericParameter("drag", "Drag coefficient", &m_dragCoefficient);
	setGroup(DRAG_COEFFICIENT, "Fluid Model|Drag force");
	setDescription(DRAG_COEFFICIENT, "Coefficient for the drag force computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(DRAG_COEFFICIENT));
	rparam->setMinValue(0.0);
}

void DragForce_Macklin2014::step()
{
	const Real density0 = m_model->getDensity0();

	const unsigned int numParticles = m_model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &ai = m_model->getAcceleration(i);
			const Vector3r &vi = m_model->getVelocity(i);
			ai -= m_dragCoefficient * static_cast<Real>(1.0) / m_model->getMass(i) * vi * (1.0 - m_model->getDensity(i) / density0);
		}
	}
}


void DragForce_Macklin2014::reset()
{
}

