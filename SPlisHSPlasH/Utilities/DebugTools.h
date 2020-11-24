#ifndef __DebugTools_h__
#define __DebugTools_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ParameterObject.h"

namespace SPH
{
	class DebugTools : public GenParam::ParameterObject
	{
	protected:
		bool m_determineThreadIds;
		std::vector<std::vector<unsigned int>> m_threadIds;
		bool m_determineNumNeighbors;
		std::vector<std::vector<unsigned int>> m_numNeighbors;
		bool m_determineVelocityChanges;
		std::vector<std::vector<Vector3r>> m_vOld;
		std::vector<std::vector<Vector3r>> m_velocityChanges;

		virtual void initParameters();

		void determineThreadIds();
		void determineNumNeighbors();
		void determineVelocityChanges();

	public:
		static int DETERMINE_THREAD_IDS;
		static int DETERMINE_NUM_NEIGHBORS;
		static int DETERMINE_VELOCITY_CHANGES;

		DebugTools();
		~DebugTools();

		void init();
		void cleanup();

		void step();
		void reset();

		void performNeighborhoodSearchSort();
		void emittedParticles(FluidModel* model, const unsigned int startIndex);
	};
}

#endif

