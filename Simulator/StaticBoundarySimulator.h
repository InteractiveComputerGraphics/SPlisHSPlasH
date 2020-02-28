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
	};
}
 
#endif
