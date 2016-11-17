#ifndef __ViscosityBase_h__
#define __ViscosityBase_h__

#include "Common.h"
#include "FluidModel.h"

namespace SPH
{
	/** \brief Base class for all viscosity methods.
	*/
	class ViscosityBase
	{
	protected:
		FluidModel *m_model;

	public:
		ViscosityBase(FluidModel *model);
		virtual ~ViscosityBase(void);

		virtual void step() = 0;
		virtual void reset() {};

		FluidModel *getModel() { return m_model; }
	};
}

#endif
