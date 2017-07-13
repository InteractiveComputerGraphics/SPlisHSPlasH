#ifndef __SimulationDataPCISPH_h__
#define __SimulationDataPCISPH_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Predictive-corrective Incompressible SPH introduced
	* by Solenthaler and Pajarola \cite Solenthaler:2009.
	*/
	class SimulationDataPCISPH
	{
		public:
			SimulationDataPCISPH();
			virtual ~SimulationDataPCISPH();

		protected:	
			FluidModel *m_model;
			Real m_pcisph_factor;

			std::vector<Vector3r> m_lastX;
			std::vector<Vector3r> m_lastV;
			std::vector<Real> m_densityAdv;
			std::vector<Real> m_pressure;
			std::vector<Vector3r> m_pressureAccel;

		public:
			/** Initialize the arrays containing the particle data.
			*/
			virtual void init(FluidModel *model);

			/** Release the arrays containing the particle data.
			*/
			virtual void cleanup();

			/** Reset the particle data.
			*/
			virtual void reset();

			/** Important: First call m_model->performNeighborhoodSearchSort() 
			 * to call the z_sort of the neighborhood search.
			 */
			void performNeighborhoodSearchSort();

			Real getPCISPH_ScalingFactor() { return m_pcisph_factor; }

			void emittedParticles(const unsigned int startIndex);

			FORCE_INLINE Vector3r &getLastPosition(const unsigned int i)
			{
				return m_lastX[i];
			}

			FORCE_INLINE const Vector3r &getLastPosition(const unsigned int i) const
			{
				return m_lastX[i];
			}

			FORCE_INLINE void setLastPosition(const unsigned int i, const Vector3r &pos)
			{
				m_lastX[i] = pos;
			}

			FORCE_INLINE Vector3r &getLastVelocity(const unsigned int i)
			{
				return m_lastV[i];
			}

			FORCE_INLINE const Vector3r &getLastVelocity(const unsigned int i) const
			{
				return m_lastV[i];
			}

			FORCE_INLINE void setLastVelocity(const unsigned int i, const Vector3r &vel)
			{
				m_lastV[i] = vel;
			}

			FORCE_INLINE const Real getDensityAdv(const unsigned int i) const
			{
				return m_densityAdv[i];
			}

			FORCE_INLINE Real& getDensityAdv(const unsigned int i)
			{
				return m_densityAdv[i];
			}

			FORCE_INLINE void setDensityAdv(const unsigned int i, const Real d)
			{
				m_densityAdv[i] = d;
			}

			FORCE_INLINE const Real getPressure(const unsigned int i) const
			{
				return m_pressure[i];
			}

			FORCE_INLINE Real& getPressure(const unsigned int i)
			{
				return m_pressure[i];
			}

			FORCE_INLINE void setPressure(const unsigned int i, const Real p)
			{
				m_pressure[i] = p;
			}

			FORCE_INLINE Vector3r &getPressureAccel(const unsigned int i)
			{
				return m_pressureAccel[i];
			}

			FORCE_INLINE const Vector3r &getPressureAccel(const unsigned int i) const
			{
				return m_pressureAccel[i];
			}

			FORCE_INLINE void setPressureAccel(const unsigned int i, const Vector3r &val)
			{
				m_pressureAccel[i] = val;
			}

	};
}

#endif