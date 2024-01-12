#ifndef __TimeStepWCSPH_h__
#define __TimeStepWCSPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataWCSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataWCSPH;

	/** \brief This class implements the Weakly Compressible SPH for Free Surface Flows approach introduced
	* by Becker and Teschner [BT07].
	*
	* References: 
	* - [BT07] Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In ACM SIGGRAPH/Eurographics Symposium on Computer Animation, SCA '07, 209-217. Aire-la-Ville, Switzerland, Switzerland, 2007. Eurographics Association. URL: http://dl.acm.org/citation.cfm?id=1272690.1272719
	*/
	class TimeStepWCSPH : public TimeStep
	{
	protected:
		Real m_stiffness;
		Real m_exponent;

		SimulationDataWCSPH m_simulationData;

		/** Determine the pressure accelerations when the pressure is already known. */
		void computePressureAccels(const unsigned int fluidModelIndex);

		virtual void performNeighborhoodSearchSort();

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);
		virtual void initParameters();

	public:
		static int STIFFNESS;
		static int EXPONENT;

		TimeStepWCSPH();
		virtual ~TimeStepWCSPH(void);

		virtual void step();
		virtual void reset();
		virtual void resize();
	};
}

#endif
