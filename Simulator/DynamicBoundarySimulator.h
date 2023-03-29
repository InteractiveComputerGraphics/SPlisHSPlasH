#ifndef __DynamicBoundarySimulator_h__
#define __DynamicBoundarySimualtor_h__
#include "BoundarySimulator.h"
namespace SPH {
	class SimulatorBase;
	class TriangleMesh;

	// This class is used for the strong coupling method. See DynamicRigidBody.
	class DynamicBoundarySimulator : public BoundarySimulator {
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
	};

}

#endif
