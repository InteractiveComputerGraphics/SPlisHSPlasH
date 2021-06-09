#ifndef __StaticBoundarySimulator_h__
#define __StaticBoundarySimulator_h__

#include "BoundarySimulator.h"

namespace SPH
{
	class SimulatorBase;
	class TriangleMesh;

	class StaticBoundarySimulator : public BoundarySimulator
	{
	protected:
		SimulatorBase *m_base;

		void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale);

	public:
		StaticBoundarySimulator(SimulatorBase *base);
		virtual ~StaticBoundarySimulator();

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
