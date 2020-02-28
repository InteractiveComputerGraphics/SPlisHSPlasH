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
		virtual void timeStep();
		virtual void initBoundaryData();
		virtual void reset();

		PBDWrapper *getPBDWrapper() { return m_pbdWrapper; }
	};
}
 
#endif
