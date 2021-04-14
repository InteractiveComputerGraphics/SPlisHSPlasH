#ifndef __PBDBoundarySimulator_h__
#define __PBDBoundarySimulator_h__

#include "Simulator/BoundarySimulator.h"
#include "PBDWrapper.h"
#include "Simulator/SimulatorBase.h"

namespace SPH
{
	class PBDBoundarySimulator : public BoundarySimulator
	{
	protected:
		PBDWrapper *m_pbdWrapper;
		SimulatorBase *m_base;

	public:
		PBDBoundarySimulator(SimulatorBase *base);
		virtual ~PBDBoundarySimulator();

		virtual void init();
		/** This function is called after the simulation scene is loaded and all
		* parameters are initialized. While reading a scene file several parameters
		* can change. The deferred init function should initialize all values which
		* depend on these parameters.
		*/
		virtual void deferredInit();
		virtual void timeStep();
		virtual void initBoundaryData();
		virtual void reset();

		PBDWrapper *getPBDWrapper() { return m_pbdWrapper; }
	};
}
 
#endif
