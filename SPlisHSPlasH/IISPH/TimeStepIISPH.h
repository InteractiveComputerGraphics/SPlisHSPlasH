#ifndef __TimeStepIISPH_h__
#define __TimeStepIISPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataIISPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataIISPH;

	/** \brief This class implements the Implicit Incompressible SPH approach introduced
	 * by Ihmsen et al. \cite Ihmsen:2014.
	 */
	class TimeStepIISPH : public TimeStep
	{
	protected:
		SimulationDataIISPH m_simulationData;
		unsigned int m_counter;

		void predictAdvection();
		void pressureSolve();
		void integration();

		/** Determine the pressure accelerations when the pressure is already known. */
		void computePressureAccels();

		/** Perform the neighborhood search for all fluid particles.
		*/
		virtual void performNeighborhoodSearch();

		virtual void emittedParticles(const unsigned int startIndex);

	public:
		TimeStepIISPH(FluidModel *model);
		virtual ~TimeStepIISPH(void);

		virtual void step();
		virtual void reset();

		const SimulationDataIISPH &getSimulationData() { return m_simulationData; };
	};
}

#endif
