#ifndef __TimeStepPCISPH_h__
#define __TimeStepPCISPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataPCISPH.h"

namespace SPH
{
	class SimulationDataPCISPH;

	/** \brief This class implements the Predictive-corrective Incompressible SPH approach introduced
	* by Solenthaler and Pajarola \cite Solenthaler:2009.
	*/
	class TimeStepPCISPH : public TimeStep
	{
	protected:
		SimulationDataPCISPH m_simulationData;
		unsigned int m_counter;

		void pressureSolve();

		/** Perform the neighborhood search for all fluid particles.
		*/
		virtual void performNeighborhoodSearch();

		virtual void emittedParticles(const unsigned int startIndex);

	public:
		TimeStepPCISPH(FluidModel *model);
		virtual ~TimeStepPCISPH(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
