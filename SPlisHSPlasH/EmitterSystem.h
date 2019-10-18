#ifndef __EmitterSystem_h__
#define __EmitterSystem_h__

#include "Common.h"
#include <vector>
#include "Emitter.h"

namespace SPH 
{	
	class TimeStep;
	class FluidModel;

	class EmitterSystem
	{
		public:
			EmitterSystem(FluidModel *model);
			virtual ~EmitterSystem();

		protected:
			FluidModel *m_model;
			static const unsigned int m_maxParticlesToReusePerStep = 50000;
			bool m_reuseParticles;
			Vector3r m_boxMin;
			Vector3r m_boxMax;
			unsigned int m_numberOfEmittedParticles;
			unsigned int m_numReusedParticles;
			std::vector <unsigned int> m_reusedParticles;
			std::vector<Emitter*> m_emitters;

			void reuseParticles();

		public:
			void enableReuseParticles(const Vector3r &boxMin = Vector3r(-1, -1, -1), const Vector3r &boxMax = Vector3r(1, 1, 1));
			void disableReuseParticles();
			void addEmitter(const unsigned int width, const unsigned int height,
				const Vector3r &pos, const Matrix3r & rotation,
				const Real velocity,
				const unsigned int type);
			unsigned int numEmitters() const { return static_cast<unsigned int>(m_emitters.size()); }
			std::vector<Emitter*> &getEmitters() { return m_emitters; }

			unsigned int numReusedParticles() const { return m_numReusedParticles; }
			unsigned int numEmittedParticles() const { return m_numberOfEmittedParticles; }

			void step();
			void reset();

			void saveState(BinaryFileWriter &binWriter);
			void loadState(BinaryFileReader &binReader);
	};
}

#endif