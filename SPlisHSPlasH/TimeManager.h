#ifndef __TimeManager_h__
#define __TimeManager_h__

#include "Common.h"
#include "Utilities/BinaryFileReaderWriter.h"
#include "ParameterObject.h"

namespace SPH
{
	/** \brief Class to manage the current simulation time and the time step size. 
	* This class is a singleton.
	*/
	class TimeManager : public GenParam::ParameterObject
	{
	private:
		Real time;
		static TimeManager *current;
		Real h;

	public:
		static int TIME_STEP_SIZE;

		TimeManager ();
		TimeManager(const TimeManager&) = delete;
		TimeManager& operator=(const TimeManager&) = delete;
		~TimeManager ();

		virtual void initParameters();

		// Singleton
		static TimeManager* getCurrent ();
		static void setCurrent (TimeManager* tm);
		static bool hasCurrent();

		Real getTime();
		void setTime(Real t);
		Real getTimeStepSize();
		void setTimeStepSize(Real tss);

		void saveState(BinaryFileWriter &binWriter);
		void loadState(BinaryFileReader &binReader);
	};
}

#endif
