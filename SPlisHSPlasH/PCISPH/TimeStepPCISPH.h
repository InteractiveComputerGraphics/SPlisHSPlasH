#ifndef __TimeStepPCISPH_h__
#define __TimeStepPCISPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataPCISPH.h"

namespace SPH
{
	class SimulationDataPCISPH;

	/** \brief This class implements the Predictive-corrective Incompressible SPH approach introduced
	* by Solenthaler and Pajarola [SP09].
	*
	* References:
	* - [SP09] B. Solenthaler and R. Pajarola. Predictive-corrective incompressible SPH. ACM Trans. Graph., 28(3):40:1-40:6, July 2009. URL: http://doi.acm.org/10.1145/1531326.1531346
	*/
	class TimeStepPCISPH : public TimeStep
	{
	protected:
		SimulationDataPCISPH m_simulationData;
		unsigned int m_iterations;
		Real m_maxError;
		unsigned int m_minIterations;
		unsigned int m_maxIterations;

		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);

		virtual void performNeighborhoodSearchSort();

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		/** Init all generic parameters */
		virtual void initParameters();

	public:
		static std::string METHOD_NAME;
		static int SOLVER_ITERATIONS;
		static int MIN_ITERATIONS;
		static int MAX_ITERATIONS;
		static int MAX_ERROR;

		TimeStepPCISPH();
		virtual ~TimeStepPCISPH(void);

		virtual void step();
		virtual void reset();
		virtual void resize();
		virtual std::string getMethodName() { return METHOD_NAME; }
		virtual int getNumIterations() { return m_iterations; }
	};
}

#endif
