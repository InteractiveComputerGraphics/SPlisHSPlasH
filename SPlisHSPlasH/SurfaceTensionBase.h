#ifndef __SurfaceTensionBase_h__
#define __SurfaceTensionBase_h__

#include "Common.h"
#include "FluidModel.h"

namespace SPH
{
	/** \brief Base class for all surface tension methods. 
	*/
	class SurfaceTensionBase
	{
	protected:
		FluidModel *m_model;

	public:
		SurfaceTensionBase(FluidModel *model);
		virtual ~SurfaceTensionBase(void);

		virtual void step() = 0;
		virtual void reset() {};

		virtual void performNeighborhoodSearchSort() {};

		FluidModel *getModel() { return m_model; }
	};
}

#endif
