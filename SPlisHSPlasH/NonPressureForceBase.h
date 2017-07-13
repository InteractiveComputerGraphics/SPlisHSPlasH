#ifndef __NonPressureForceBaseBase_h__
#define __NonPressureForceBaseBase_h__

#include "Common.h"
#include "FluidModel.h"

namespace SPH
{
	/** \brief Base class for all non-pressure force methods.
	*/
	class NonPressureForceBase
	{
	protected:
		FluidModel *m_model;

	public:
		NonPressureForceBase(FluidModel *model);
		virtual ~NonPressureForceBase(void);

		virtual void step() = 0;
		virtual void reset() {};

		virtual void performNeighborhoodSearchSort() {};
		virtual void emittedParticles(const unsigned int startIndex) {};

		FluidModel *getModel() { return m_model; }
	};
}

#endif
