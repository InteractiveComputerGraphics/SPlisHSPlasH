#ifndef __TimeStepDFSPH_h__
#define __TimeStepDFSPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataDFSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

#define USE_WARMSTART
#define USE_WARMSTART_V

namespace SPH
{
	class SimulationDataDFSPH;

	/** \brief This class implements the Divergence-free Smoothed Particle Hydrodynamics approach introduced
	* by Bender and Koschier \cite Bender:2015, \cite Bender2017, \cite KBST19.
	*/
	class TimeStepDFSPH : public TimeStep
	{
	protected:
		SimulationDataDFSPH m_simulationData;
		unsigned int m_counter;
		const Real m_eps = 1.0e-5;
		bool m_enableDivergenceSolver;
		unsigned int m_iterationsV;
		Real m_maxErrorV;
		unsigned int m_maxIterationsV;

		void computeDFSPHFactor(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void divergenceSolve();
		void divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int index, const int numParticles, const Real h, const Real density0);
		void computeDensityChange(const unsigned int fluidModelIndex, const unsigned int index, const Real h);

#ifdef USE_WARMSTART_V
		void warmstartDivergenceSolve(const unsigned int fluidModelIndex);
#endif
#ifdef USE_WARMSTART
		void warmstartPressureSolve(const unsigned int fluidModelIndex);
#endif

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();
		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		virtual void initParameters();

	public:
		static int SOLVER_ITERATIONS_V;
		static int MAX_ITERATIONS_V;
		static int MAX_ERROR_V;
		static int USE_DIVERGENCE_SOLVER;

		TimeStepDFSPH();
		virtual ~TimeStepDFSPH(void);

		virtual void step();
		virtual void reset();

		virtual void resize();
	};
}

#endif
