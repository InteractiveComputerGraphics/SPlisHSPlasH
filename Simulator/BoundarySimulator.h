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
		virtual void timeStep() {}
		virtual void initBoundaryData() {}
		virtual void reset() {}

		void updateBoundaryForces();		
	};
}
 
#endif
