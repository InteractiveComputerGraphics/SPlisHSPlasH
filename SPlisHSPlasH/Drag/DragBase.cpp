#include "DragBase.h"

using namespace SPH;

DragBase::DragBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_dragCoefficient = 0.01;
}

DragBase::~DragBase(void)
{
}



