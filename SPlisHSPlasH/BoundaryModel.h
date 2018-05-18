#ifndef __BoundaryModel_h__
#define __BoundaryModel_h__

#include "Common.h"
#include <vector>

#include "RigidBodyObject.h"
#include "SPHKernels.h"

namespace SPH 
{	
	class TimeStep;

	/** \brief The fluid model stores the particle and simulation information 
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
			std::vector<Real> m_boundaryPsi;
			std::vector<Vector3r> m_f;
			bool m_sorted;
			unsigned int m_pointSetIndex;

		public:
			unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }

			void computeBoundaryPsi(const Real density0);

			virtual void reset();

			void performNeighborhoodSearchSort();

			void initModel(RigidBodyObject *rbo, const unsigned int numBoundaryParticles, Vector3r *boundaryParticles);
			RigidBodyObject* getRigidBodyObject() { return m_rigidBody; }
			
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
			
			FORCE_INLINE Vector3r &getForce(const unsigned int i)
			{
				return m_f[i];
			}

			FORCE_INLINE const Vector3r &getForce(const unsigned int i) const
			{
				return m_f[i];
			}

			FORCE_INLINE void setForce(const unsigned int i, const Vector3r &f)
			{
				m_f[i] = f;
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

			FORCE_INLINE const Real& getBoundaryPsi(const unsigned int i) const
			{
				return m_boundaryPsi[i];
			}

			FORCE_INLINE Real& getBoundaryPsi(const unsigned int i)
			{
				return m_boundaryPsi[i];
			}

			FORCE_INLINE void setBoundaryPsi(const unsigned int i, const Real &val)
			{
				m_boundaryPsi[i] = val;
			}
	};
}

#endif