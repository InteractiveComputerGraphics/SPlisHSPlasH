#ifndef __DragForce_Macklin2014_h__
#define __DragForce_Macklin2014_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "DragBase.h"

namespace SPH
{
	/** \brief This class implements the drag force computation introduced
	* by Macklin et al. \cite Macklin:2014.
	*/
	class DragForce_Macklin2014 : public DragBase
	{
	public:
		DragForce_Macklin2014(FluidModel *model);
		virtual ~DragForce_Macklin2014(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
