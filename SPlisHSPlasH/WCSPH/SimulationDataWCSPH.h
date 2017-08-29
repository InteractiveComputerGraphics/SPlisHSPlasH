#ifndef __SimulationDataWCSPH_h__
#define __SimulationDataWCSPH_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Weakly Compressible SPH for Free Surface Flows introduced
	* by Becker and Teschner \cite Becker:2007.
	*/
	class SimulationDataWCSPH
	{
		public:
			SimulationDataWCSPH();
			virtual ~SimulationDataWCSPH();

		protected:	
			FluidModel *m_model;

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

			void emittedParticles(const unsigned int startIndex);

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