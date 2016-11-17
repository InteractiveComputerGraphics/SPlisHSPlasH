#ifndef __TimeStepDFSPH_h__
#define __TimeStepDFSPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataDFSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataDFSPH;

	/** \brief This class implements the Divergence-free Smoothed Particle Hydrodynamics approach introduced
	* by Bender and Koschier \cite Bender:2015, \cite Bender2016.
	*/
	class TimeStepDFSPH : public TimeStep
	{
	protected:
		SimulationDataDFSPH m_simulationData;
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
		TimeStepDFSPH(FluidModel *model);
		virtual ~TimeStepDFSPH(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
