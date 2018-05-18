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

		void predictAdvection(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void integration(const unsigned int fluidModelIndex);

		/** Determine the pressure accelerations when the pressure is already known. */
		void computePressureAccels(const unsigned int fluidModelIndex);

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

	public:
		TimeStepIISPH();
		virtual ~TimeStepIISPH(void);

		virtual void step();
		virtual void reset();
		virtual void resize();

		const SimulationDataIISPH &getSimulationData() { return m_simulationData; };
	};
}

#endif
