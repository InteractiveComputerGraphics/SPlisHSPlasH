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
	 * by Ihmsen et al. [ICS+14].
	*
	* References:
	* - [ICS+14] Markus Ihmsen, Jens Cornelis, Barbara Solenthaler, Christopher Horvath, and Matthias Teschner. Implicit incompressible SPH. IEEE Transactions on Visualization and Computer Graphics, 20(3):426-435, March 2014. URL: http://dx.doi.org/10.1109/TVCG.2013.105
	*/
	class TimeStepIISPH : public TimeStep
	{
	protected:
		SimulationDataIISPH m_simulationData;
		unsigned int m_iterations;
		Real m_maxError;
		unsigned int m_minIterations;
		unsigned int m_maxIterations;

		void predictAdvection(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void integration(const unsigned int fluidModelIndex);

		/** Determine the pressure accelerations when the pressure is already known. */
		void computePressureAccels(const unsigned int fluidModelIndex);

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

		TimeStepIISPH();
		virtual ~TimeStepIISPH(void);

		virtual void step();
		virtual void reset();
		virtual void resize();
		virtual std::string getMethodName() { return METHOD_NAME; }

		const SimulationDataIISPH &getSimulationData() { return m_simulationData; };
		virtual int getNumIterations() { return m_iterations; }
	};
}

#endif
