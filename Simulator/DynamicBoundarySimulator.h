#ifndef __DynamicBoundarySimulator_h__
#define __DynamicBoundarySimualtor_h__
#include "BoundarySimulator.h"
namespace SPH {
	class SimulatorBase;
	class TriangleMesh;

	class DynamicBoundarySimulator : public BoundarySimulator {
	private:
		Real m_dampingCoeff = 0.0;


	protected:
		SimulatorBase* m_base;

	public:
		DynamicBoundarySimulator(SimulatorBase* base);
		virtual ~DynamicBoundarySimulator();

		virtual void initBoundaryData();
		/** This function is called after the simulation scene is loaded and all
		* parameters are initialized. While reading a scene file several parameters
		* can change. The deferred init function should initialize all values which
		* depend on these parameters.
		*/
		virtual void deferredInit();

		virtual void timeStep();
		virtual void reset();

		FORCE_INLINE const Real& getDampingCoeff() const {
			return m_dampingCoeff;
		}
		FORCE_INLINE Real& getDampingCoeff() {
			return m_dampingCoeff;
		}
		FORCE_INLINE void setDampingCoeff(const Real& value) {
			m_dampingCoeff = value;
		}

	};

}

#endif
