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
	* by Macklin and Mueller [MM13,BMO+14,BMM15].
	*
	* References:
	* - [MM13] Miles Macklin and Matthias Müller. Position based fluids. ACM Trans. Graph., 32(4):104:1-104:12, July 2013. URL: http://doi.acm.org/10.1145/2461912.2461984
	* - [BMO+14] Jan Bender, Matthias Müller, Miguel A. Otaduy, Matthias Teschner, and Miles Macklin. A survey on position-based simulation methods in computer graphics. Computer Graphics Forum, 33(6):228-251, 2014. URL: http://dx.doi.org/10.1111/cgf.12346
	* - [BMM15] Jan Bender, Matthias Müller, and Miles Macklin. Position-based simulation methods in computer graphics. In EUROGRAPHICS 2015 Tutorials. Eurographics Association, 2015. URL: http://dx.doi.org/10.2312/egt.20151045
	*/
	class TimeStepPBF : public TimeStep
	{
	protected:
		SimulationDataPBF m_simulationData;
		int m_velocityUpdateMethod;

		/** Perform a position-based correction step for the following density constraint:\n
		*  \f$C(\mathbf{x}) = \left (\frac{\rho_i}{\rho_0} - 1 \right )= 0\f$\n
		*/
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);

		virtual void performNeighborhoodSearchSort();

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
