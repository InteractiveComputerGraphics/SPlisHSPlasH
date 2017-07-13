#ifndef __Emitter_h__
#define __Emitter_h__

#include "Common.h"
#include <vector>


namespace SPH 
{	
	class TimeStep;

	class Emitter 
	{
		public:
			Emitter(const unsigned int width, const unsigned int height,
				const Vector3r &pos, const Vector3r &dir, const Vector3r &initialVel,
				const Real emitsPerSecond, const unsigned int type = 0);
			virtual ~Emitter();

	protected:
			unsigned int m_width; 
			unsigned int m_height;
			Vector3r m_x;
			Vector3r m_dir;
			Vector3r m_v;
			Real m_emitsPerSecond;
			unsigned int m_type;
			Real m_nextEmitTime;
			unsigned int m_emitCounter;

		public:
			void emitParticles(TimeStep *timeStep, std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles);
			void emitParticlesCircle(TimeStep *timeStep, std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles);
			Real getNextEmitTime() const { return m_nextEmitTime; }
			void setNextEmitTime(Real val) { m_nextEmitTime = val; }
			static void getOrthogonalVectors(const Vector3r &vec, Vector3r &x, Vector3r &y);

			void step(TimeStep *timeStep, std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles);
			virtual void reset();

	};
}

#endif