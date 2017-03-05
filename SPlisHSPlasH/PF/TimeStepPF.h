#ifndef __TimeStepPF_h__
#define __TimeStepPF_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataPF.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	/** \brief This class implements the Projective Fluids approach introduced
	* by Weiler, Koschier and Bender \cite Weiler:2016.
	*/
	class TimeStepPF : public TimeStep
	{
	protected:
		SimulationDataPF m_simulationData;
		unsigned int m_counter;

		void initialGuessForPositions();
		void countFluidNeighbors();
		void solvePDConstraints();
		void updateVelocity();

		/** Perform the neighborhood search for all fluid particles.
		*/
		virtual void performNeighborhoodSearch();

	public:
		TimeStepPF(FluidModel *model);
		virtual ~TimeStepPF(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
