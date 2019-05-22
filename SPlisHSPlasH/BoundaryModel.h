#ifndef __BoundaryModel_h__
#define __BoundaryModel_h__

#include "Common.h"
#include <vector>

#include "RigidBodyObject.h"
#include "SPHKernels.h"

namespace SPH 
{	
	class TimeStep;

	/** \brief The boundary model stores the information required for boundary handling
	*/
	class BoundaryModel 
	{
		public:
			BoundaryModel();
			virtual ~BoundaryModel();

		protected:
			RigidBodyObject *m_rigidBody;
			std::vector<Vector3r> m_x0;
			std::vector<Vector3r> m_x;
			std::vector<Vector3r> m_v;
			std::vector<Real> m_V;
			std::vector<Vector3r> m_forcePerThread;
			std::vector<Vector3r> m_torquePerThread;
			bool m_sorted;
			unsigned int m_pointSetIndex;

		public:
			unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }

			void computeBoundaryVolume();

			virtual void reset();

			void performNeighborhoodSearchSort();

			void initModel(RigidBodyObject *rbo, const unsigned int numBoundaryParticles, Vector3r *boundaryParticles);
			RigidBodyObject* getRigidBodyObject() { return m_rigidBody; }

			FORCE_INLINE void addForce(const Vector3r &pos, const Vector3r &f)
			{
				if (m_rigidBody->isDynamic())
				{
					#ifdef _OPENMP
					int tid = omp_get_thread_num();
					#else
					int tid = 0;
					#endif
					m_forcePerThread[tid] += f;
					m_torquePerThread[tid] += (pos - m_rigidBody->getPosition()).cross(f);
				}
			}

			void getForceAndTorque(Vector3r &force, Vector3r &torque);
			void clearForceAndTorque();

			FORCE_INLINE Vector3r &getPosition0(const unsigned int i)
			{
				return m_x0[i];
			}

			FORCE_INLINE const Vector3r &getPosition0(const unsigned int i) const
			{
				return m_x0[i];
			}

			FORCE_INLINE void setPosition0(const unsigned int i, const Vector3r &pos)
			{
				m_x0[i] = pos;
			}

			FORCE_INLINE Vector3r &getPosition(const unsigned int i)
			{
				return m_x[i];
			}

			FORCE_INLINE const Vector3r &getPosition(const unsigned int i) const
			{
				return m_x[i];
			}

			FORCE_INLINE void setPosition(const unsigned int i, const Vector3r &pos)
			{
				m_x[i] = pos;
			}

			FORCE_INLINE Vector3r &getVelocity(const unsigned int i)
			{
				return m_v[i];
			}

			FORCE_INLINE const Vector3r &getVelocity(const unsigned int i) const
			{
				return m_v[i];
			}

			FORCE_INLINE void setVelocity(const unsigned int i, const Vector3r &vel)
			{
				m_v[i] = vel;
			}

			FORCE_INLINE const Real& getVolume(const unsigned int i) const
			{
				return m_V[i];
			}

			FORCE_INLINE Real& getVolume(const unsigned int i)
			{
				return m_V[i];
			}

			FORCE_INLINE void setVolume(const unsigned int i, const Real &val)
			{
				m_V[i] = val;
			}
	};
}

#endif