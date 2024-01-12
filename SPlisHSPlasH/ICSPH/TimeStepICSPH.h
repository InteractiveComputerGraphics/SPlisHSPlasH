#ifndef __TimeStepICSPH_h__
#define __TimeStepICSPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataICSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataICSPH;

	/** \brief This class implements the Implicit Compressible SPH approach introduced
	 * by Gissler et al. [GHB+20].
	*
	* References:
	* - [GHB+20] Christoph Gissler, Andreas Henne, Stefan Band, Andreas Peer and Matthias Teschner. An Implicit Compressible SPH Solver for Snow Simulation. ACM Transactions on Graphics, 39(4). URL: https://doi.org/10.1145/3386569.3392431
	*/
	class TimeStepICSPH : public TimeStep
	{
	protected:
		SimulationDataICSPH m_simulationData;
		Real m_lambda;
		bool m_clamping;
		const Real m_psi = 1.5;

		void computeDensityAdv(const unsigned int fluidModelIndex);
		void compute_aii(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void integration(const unsigned int fluidModelIndex);

		/** Determine the pressure accelerations when the pressure is already known. */
		void computePressureAccels(const unsigned int fluidModelIndex);

		virtual void performNeighborhoodSearchSort();

		virtual void initParameters();
		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

	public:
		static int LAMBDA;
		static int PRESSURE_CLAMPING;

		TimeStepICSPH();
		virtual ~TimeStepICSPH(void);

		virtual void step();
		virtual void reset();
		virtual void resize();

		const SimulationDataICSPH &getSimulationData() { return m_simulationData; };
	};
}

#endif
