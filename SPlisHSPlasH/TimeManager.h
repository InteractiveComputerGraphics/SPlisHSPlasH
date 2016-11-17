#ifndef __TimeManager_h__
#define __TimeManager_h__

#include "Common.h"

namespace SPH
{
	/** \brief Class to manage the current simulation time and the time step size. 
	* This class is a singleton.
	*/
	class TimeManager
	{
	private:
		Real time;
		static TimeManager *current;
		Real h;

	public:
		TimeManager ();
		~TimeManager ();

		// Singleton
		static TimeManager* getCurrent ();
		static void setCurrent (TimeManager* tm);
		static bool hasCurrent();

		Real getTime();
		void setTime(Real t);
		Real getTimeStepSize();
		void setTimeStepSize(Real tss);
	};
}

#endif
