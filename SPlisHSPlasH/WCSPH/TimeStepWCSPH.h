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
	* by Becker and Teschner \cite Becker:2007.
	*/
	class TimeStepWCSPH : public TimeStep
	{
	protected:
		Real m_stiffness;
		Real m_exponent;

		SimulationDataWCSPH m_simulationData;
		unsigned int m_counter;

		/** Determine the pressure accelerations when the pressure is already known. */
		void computePressureAccels();

		/** Perform the neighborhood search for all fluid particles.
		*/
		virtual void performNeighborhoodSearch();

		virtual void emittedParticles(const unsigned int startIndex);

	public:
		TimeStepWCSPH(FluidModel *model);
		virtual ~TimeStepWCSPH(void);

		virtual void step();
		virtual void reset();

		Real getStiffness() const { return m_stiffness; }
		void setStiffness(Real val) { m_stiffness = val; }

		Real getExponent() const { return m_exponent; }
		void setExponent(Real val) { m_exponent = val; }
	};
}

#endif
