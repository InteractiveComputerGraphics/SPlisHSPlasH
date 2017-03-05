#ifndef __TimeStepPF_h__
#define __TimeStepPF_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataPF.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataPF;

	/** \brief This class implements the Divergence-free Smoothed Particle Hydrodynamics approach introduced
	* by Bender and Koschier \cite Bender:2015, \cite Bender2016.
	*/
	class TimeStepPF : public TimeStep
	{
	protected:
		SimulationDataPF m_simulationData;
		unsigned int m_counter;

		void computeDFSPHFactor();
		void pressureSolve();
		void divergenceSolve();
		void computeDensityAdv(const unsigned int index, const int numParticles, const Real h, const Real density0);
		void computeDensityChange(const unsigned int index, const Real h, const Real density0);

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
