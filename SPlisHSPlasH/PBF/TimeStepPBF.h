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
		unsigned int m_velocityUpdateMethod;

		/** Perform a position-based correction step for the following density constraint:\n
		*  \f$C(\mathbf{x}) = \left (\frac{\rho_i}{\rho_0} - 1 \right )= 0\f$\n
		*/
		void pressureSolve();

		/** Perform the neighborhood search for all fluid particles. 
		*/
		virtual void performNeighborhoodSearch();

		virtual void emittedParticles(const unsigned int startIndex);

	public:
		/** Initialize the simulation data required for this method. */
		TimeStepPBF(FluidModel *model);
		virtual ~TimeStepPBF(void);

		/** Perform a simulation step. */
		virtual void step();

		/** Reset the simulation method. */
		virtual void reset();

		unsigned int getVelocityUpdateMethod() const { return m_velocityUpdateMethod; }
		void setVelocityUpdateMethod(unsigned int val) { m_velocityUpdateMethod = val; }
	};
}

#endif
