#include "DragBase.h"

using namespace SPH;
using namespace GenParam;

int DragBase::DRAG_COEFFICIENT = -1;


DragBase::DragBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_dragCoefficient = 0.01;
}

DragBase::~DragBase(void)
{
}

void DragBase::initParameters()
{
	NonPressureForceBase::initParameters();

	DRAG_COEFFICIENT = createNumericParameter("drag", "Drag coefficient", &m_dragCoefficient);
	setGroup(DRAG_COEFFICIENT, "Drag force");
	setDescription(DRAG_COEFFICIENT, "Coefficient for the drag force computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(DRAG_COEFFICIENT));
	rparam->setMinValue(0.0);
}


