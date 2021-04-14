#ifndef __NonPressureForceBaseBase_h__
#define __NonPressureForceBaseBase_h__

#include "Common.h"
#include "FluidModel.h"
#include "ParameterObject.h"

namespace SPH
{
	/** \brief Base class for all non-pressure force methods.
	*/
	class NonPressureForceBase : public GenParam::ParameterObject
	{
	protected:
		FluidModel *m_model;

	public:
		NonPressureForceBase(FluidModel *model);
		NonPressureForceBase(const NonPressureForceBase&) = delete;
        NonPressureForceBase& operator=(const NonPressureForceBase&) = delete;
		virtual ~NonPressureForceBase(void);

		virtual void step() = 0;
		virtual void reset() {};

		virtual void performNeighborhoodSearchSort() {};
		virtual void emittedParticles(const unsigned int startIndex) {};

		virtual void saveState(BinaryFileWriter &binWriter) {};
		virtual void loadState(BinaryFileReader &binReader) {};

		FluidModel *getModel() { return m_model; }

		virtual void init();
		/** This function is called after the simulation scene is loaded and all
		* parameters are initialized. While reading a scene file several parameters
		* can change. The deferred init function should initialize all values which
		* depend on these parameters.
		*/
		virtual void deferredInit() {};
	};
}

#endif
