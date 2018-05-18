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
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

	public:
		TimeStepPCISPH();
		virtual ~TimeStepPCISPH(void);

		virtual void step();
		virtual void reset();
		virtual void resize();
	};
}

#endif
