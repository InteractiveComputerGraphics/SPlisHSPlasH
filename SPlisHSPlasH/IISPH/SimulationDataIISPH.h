#ifndef __SimulationDataIISPH_h__
#define __SimulationDataIISPH_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Implicit Incompressible SPH introduced
	* by Ihmsen et al. \cite Ihmsen:2014.
	*/
	class SimulationDataIISPH
	{
		public:
			SimulationDataIISPH();
			virtual ~SimulationDataIISPH();

		protected:	
			FluidModel *m_model;

			std::vector<Real> m_aii;
			std::vector<Vector3r> m_dii;
			std::vector<Vector3r> m_dij_pj;
			std::vector<Real> m_density_adv;
			std::vector<Real> m_pressure;
			std::vector<Real> m_lastPressure;
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

			void emittedParticles(const unsigned int startIndex);

			FORCE_INLINE const Real getAii(const unsigned int i) const
			{
				return m_aii[i];
			}

			FORCE_INLINE Real& getAii(const unsigned int i)
			{
				return m_aii[i];
			}

			FORCE_INLINE void setAii(const unsigned int i, const Real aii)
			{
				m_aii[i] = aii;
			}

			FORCE_INLINE Vector3r &getDii(const unsigned int i)
			{
				return m_dii[i];
			}

			FORCE_INLINE const Vector3r &getDii(const unsigned int i) const
			{
				return m_dii[i];
			}

			FORCE_INLINE void setDii(const unsigned int i, const Vector3r &val)
			{
				m_dii[i] = val;
			}

			FORCE_INLINE Vector3r &getDij_pj(const unsigned int i)
			{
				return m_dij_pj[i];
			}

			FORCE_INLINE const Vector3r &getDij_pj(const unsigned int i) const
			{
				return m_dij_pj[i];
			}

			FORCE_INLINE void setDij_pj(const unsigned int i, const Vector3r &val)
			{
				m_dij_pj[i] = val;
			}

			FORCE_INLINE const Real getDensityAdv(const unsigned int i) const
			{
				return m_density_adv[i];
			}

			FORCE_INLINE Real& getDensityAdv(const unsigned int i)
			{
				return m_density_adv[i];
			}

			FORCE_INLINE void setDensityAdv(const unsigned int i, const Real d)
			{
				m_density_adv[i] = d;
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

			FORCE_INLINE const Real getLastPressure(const unsigned int i) const
			{
				return m_lastPressure[i];
			}

			FORCE_INLINE Real& getLastPressure(const unsigned int i)
			{
				return m_lastPressure[i];
			}

			FORCE_INLINE void setLastPressure(const unsigned int i, const Real p)
			{
				m_lastPressure[i] = p;
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