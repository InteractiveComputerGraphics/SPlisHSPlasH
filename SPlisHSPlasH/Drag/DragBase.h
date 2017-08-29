#ifndef __DragBase_h__
#define __DragBase_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief Base class for all drag force methods.
	*/
	class DragBase : public NonPressureForceBase
	{
	protected:
		Real m_dragCoefficient;

	public:
		DragBase(FluidModel *model);
		virtual ~DragBase(void);

		Real getDragCoefficient() const { return m_dragCoefficient; }
		void setDragCoefficient(Real val) { m_dragCoefficient = val; }
	};
}

#endif
