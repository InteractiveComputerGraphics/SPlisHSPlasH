#ifndef __TimeStepPBF_h__
#define __TimeStepPBF_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataPBF.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataPBF;


	/** \brief This class implements the position-based fluids approach introduced 
	* by Macklin and Mueller \cite Macklin:2013:PBF, \cite BMOTM2014, \cite BMM2015.
	*/
	class TimeStepPBF : public TimeStep
	{
	protected:
		SimulationDataPBF m_simulationData;
		unsigned int m_counter;
		int m_velocityUpdateMethod;

		/** Perform a position-based correction step for the following density constraint:\n
		*  \f$C(\mathbf{x}) = \left (\frac{\rho_i}{\rho_0} - 1 \right )= 0\f$\n
		*/
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);

		/** Perform the neighborhood search for all fluid particles. 
		*/
		void performNeighborhoodSearch();

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		virtual void initParameters();

	public:
		static int VELOCITY_UPDATE_METHOD;
		static int ENUM_PBF_FIRST_ORDER;
		static int ENUM_PBF_SECOND_ORDER;

		/** Initialize the simulation data required for this method. */
		TimeStepPBF();
		virtual ~TimeStepPBF(void);

		/** Perform a simulation step. */
		virtual void step();

		/** Reset the simulation method. */
		virtual void reset();
		virtual void resize();
	};
}

#endif
