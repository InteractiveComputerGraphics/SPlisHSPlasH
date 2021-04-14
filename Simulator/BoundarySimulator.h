#ifndef __BoundarySimulator_h__
#define __BoundarySimulator_h__

#include "SPlisHSPlasH/Common.h"

namespace SPH
{
	class BoundarySimulator 
	{
	public:
		BoundarySimulator() {}
		virtual ~BoundarySimulator() {}
		virtual void init() {}
		/** This function is called after the simulation scene is loaded and all
		* parameters are initialized. While reading a scene file several parameters
		* can change. The deferred init function should initialize all values which
		* depend on these parameters.
		*/
		virtual void deferredInit() {}
		virtual void timeStep() {}
		virtual void initBoundaryData() {}
		virtual void reset() {}

		void updateBoundaryForces();		
	};
}
 
#endif
