#ifndef __FluidModel_h__
#define __FluidModel_h__

#include "Common.h"
#include <vector>

#include "CompactNSearch.h"
#include "RigidBodyObject.h"
#include "SPHKernels.h"
#include "EmitterSystem.h"

namespace SPH 
{	
	class TimeStep;

	/** \brief The fluid model stores the particle and simulation information 
	*/
	class FluidModel 
	{
		public:
			FluidModel();
			virtual ~FluidModel();

			/** \brief Struct to store the state of a particle object (x0, x, v).
			*/
			struct ParticleObject
			{
				std::vector<Vector3r> m_x0;
				std::vector<Vector3r> m_x;
				std::vector<Vector3r> m_v;				
				unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }
			};

			/** \brief Struct to store the pseudo masses and forces of the sampling of a rigid body object.
			*/
			struct RigidBodyParticleObject : public ParticleObject
			{
				RigidBodyObject *m_rigidBody;
				std::vector<Real> m_boundaryPsi;
				std::vector<Vector3r> m_f;
			};

			typedef PrecomputedKernel<CubicKernel, 10000> PrecomputedCubicKernel;

	protected:
			EmitterSystem m_emitterSystem;
			Vector3r m_gravitation;
			unsigned int m_kernelMethod;
			unsigned int m_gradKernelMethod;
			Real m_W_zero;
			Real(*m_kernelFct)(const Vector3r &);
			Vector3r(*m_gradKernelFct)(const Vector3r &r);

			std::vector<ParticleObject*> m_particleObjects;

			// Mass
			// If the mass is zero, the particle is static
			std::vector<Real> m_masses;
			std::vector<Vector3r> m_a;
			std::vector<Vector3r> m_v0;

			// initial position
			std::vector<Real> m_density;

			Real m_density0;
			Real m_particleRadius;
			Real m_supportRadius;
			CompactNSearch::NeighborhoodSearch *m_neighborhoodSearch;

			unsigned int m_numActiveParticles;
			unsigned int m_numActiveParticles0;

			void initMasses();
			void computeBoundaryPsi(const unsigned int body);

			/** Resize the arrays containing the particle data.
			*/
			virtual void resizeFluidParticles(const unsigned int newSize);
			
			/** Release the arrays containing the particle data.
			*/
			virtual void releaseFluidParticles();

		public:
			void setNumActiveParticles(const unsigned int num);

			virtual void cleanupModel();
			virtual void reset();

			void updateBoundaryPsi();

			void initModel(const unsigned int nFluidParticles, Vector3r* fluidParticles, Vector3r* fluidVelocities, const unsigned int nMaxEmitterParticles);
			void addRigidBodyObject(RigidBodyObject *rbo, const unsigned int numBoundaryParticles, Vector3r *boundaryParticles);
			
			RigidBodyParticleObject *getRigidBodyParticleObject(const unsigned int index) { return (FluidModel::RigidBodyParticleObject*) m_particleObjects[index + 1]; }
			const unsigned int numParticles() const { return static_cast<unsigned int>(m_masses.size()); }
			unsigned int numActiveParticles() const;
			const unsigned int numberOfRigidBodyParticleObjects() const { return static_cast<unsigned int>(m_particleObjects.size()-1); }

			FORCE_INLINE Real getDensity0() const { return m_density0; }
			void setDensity0(const Real v);
			Real getSupportRadius() const { return m_supportRadius; }
			Real getParticleRadius() const { return m_particleRadius; }
			void setParticleRadius(Real val);

			unsigned int getNumActiveParticles0() const { return m_numActiveParticles0; }
			void setNumActiveParticles0(unsigned int val) { m_numActiveParticles0 = val; }

			SPH::EmitterSystem& getEmitterSystem() { return m_emitterSystem; }

			unsigned int getKernel() const { return m_kernelMethod; }
			void setKernel(unsigned int val);
			unsigned int getGradKernel() const { return m_gradKernelMethod; }
			void setGradKernel(unsigned int val);

			FORCE_INLINE Real W_zero() const { return m_W_zero; }
			FORCE_INLINE Real W(const Vector3r &r) const { return m_kernelFct(r); }
			FORCE_INLINE Vector3r gradW(const Vector3r &r) { return m_gradKernelFct(r); }

			const SPH::Vector3r& getGravitation() const { return m_gravitation; }
			void setGravitation(const SPH::Vector3r &val) { m_gravitation = val; }
			
			CompactNSearch::NeighborhoodSearch* getNeighborhoodSearch() { return m_neighborhoodSearch; }
			void performNeighborhoodSearchSort();

			FORCE_INLINE unsigned int numberOfPointSets() const
			{
				return static_cast<unsigned int>(m_neighborhoodSearch->n_point_sets());
			}

			FORCE_INLINE unsigned int numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int index) const
			{
				return static_cast<unsigned int>(m_neighborhoodSearch->point_set(0).n_neighbors(pointSetIndex, index));
			}

			FORCE_INLINE unsigned int getNeighbor(const unsigned int pointSetIndex, const unsigned int index, const unsigned int k) const
			{
				return m_neighborhoodSearch->point_set(0).neighbor(pointSetIndex, index, k);
			}

			FORCE_INLINE Vector3r &getPosition0(const unsigned int objectIndex, const unsigned int i)
			{
				return m_particleObjects[objectIndex]->m_x0[i];
			}

			FORCE_INLINE const Vector3r &getPosition0(const unsigned int objectIndex, const unsigned int i) const
			{
				return m_particleObjects[objectIndex]->m_x0[i];
			}

			FORCE_INLINE void setPosition0(const unsigned int objectIndex, const unsigned int i, const Vector3r &pos)
			{
				m_particleObjects[objectIndex]->m_x0[i] = pos;
			}

			FORCE_INLINE Vector3r &getPosition(const unsigned int objectIndex, const unsigned int i)
			{
				return m_particleObjects[objectIndex]->m_x[i];
			}

			FORCE_INLINE const Vector3r &getPosition(const unsigned int objectIndex, const unsigned int i) const
			{
				return m_particleObjects[objectIndex]->m_x[i];
			}

			FORCE_INLINE void setPosition(const unsigned int objectIndex, const unsigned int i, const Vector3r &pos)
			{
				m_particleObjects[objectIndex]->m_x[i] = pos;
			}

			FORCE_INLINE Vector3r &getVelocity(const unsigned int objectIndex, const unsigned int i)
			{
				return m_particleObjects[objectIndex]->m_v[i];
			}

			FORCE_INLINE const Vector3r &getVelocity(const unsigned int objectIndex, const unsigned int i) const
			{
				return m_particleObjects[objectIndex]->m_v[i];
			}

			FORCE_INLINE void setVelocity(const unsigned int objectIndex, const unsigned int i, const Vector3r &vel)
			{
				m_particleObjects[objectIndex]->m_v[i] = vel;
			}

			FORCE_INLINE Vector3r &getVelocity0(const unsigned int i)
			{
				return m_v0[i];
			}

			FORCE_INLINE const Vector3r &getVelocity0(const unsigned int i) const
			{
				return m_v0[i];
			}

			FORCE_INLINE void setVelocity0(const unsigned int i, const Vector3r &vel)
			{
				m_v0[i] = vel;
			}

			FORCE_INLINE Vector3r &getAcceleration(const unsigned int i)
			{
				return m_a[i];
			}

			FORCE_INLINE const Vector3r &getAcceleration(const unsigned int i) const
			{
				return m_a[i];
			}

			FORCE_INLINE void setAcceleration(const unsigned int i, const Vector3r &accel)
			{
				m_a[i] = accel;
			}

			FORCE_INLINE Vector3r &getForce(const unsigned int objectIndex, const unsigned int i)
			{
				return static_cast<RigidBodyParticleObject*>(m_particleObjects[objectIndex])->m_f[i];
			}

			FORCE_INLINE const Vector3r &getForce(const unsigned int objectIndex, const unsigned int i) const
			{
				return static_cast<RigidBodyParticleObject*>(m_particleObjects[objectIndex])->m_f[i];
			}

			FORCE_INLINE void setForce(const unsigned int objectIndex, const unsigned int i, const Vector3r &f)
			{
				static_cast<RigidBodyParticleObject*>(m_particleObjects[objectIndex])->m_f[i] = f;
			}

			FORCE_INLINE const Real getMass(const unsigned int i) const
			{
				return m_masses[i];
			}

			FORCE_INLINE Real& getMass(const unsigned int i)
			{
				return m_masses[i];
			}

			FORCE_INLINE void setMass(const unsigned int i, const Real mass)
			{
				m_masses[i] = mass;
			}

			FORCE_INLINE const Real& getBoundaryPsi(const unsigned int objectIndex, const unsigned int i) const
			{
				return static_cast<RigidBodyParticleObject*>(m_particleObjects[objectIndex])->m_boundaryPsi[i];
			}

			FORCE_INLINE Real& getBoundaryPsi(const unsigned int objectIndex, const unsigned int i)
			{
				return static_cast<RigidBodyParticleObject*>(m_particleObjects[objectIndex])->m_boundaryPsi[i];
			}

			FORCE_INLINE void setBoundaryPsi(const unsigned int objectIndex, const unsigned int i, const Real &val)
			{
				static_cast<RigidBodyParticleObject*>(m_particleObjects[objectIndex])->m_boundaryPsi[i] = val;
			}

			FORCE_INLINE const Real& getDensity(const unsigned int i) const
			{
				return m_density[i];
			}

			FORCE_INLINE Real& getDensity(const unsigned int i)
			{
				return m_density[i];
			}

			FORCE_INLINE void setDensity(const unsigned int i, const Real &val)
			{
				m_density[i] = val;
			}
	};
}

#endif