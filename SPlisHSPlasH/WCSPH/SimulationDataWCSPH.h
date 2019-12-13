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
			std::vector<std::vector<Real>> m_pressure;
			std::vector<std::vector<Vector3r>> m_pressureAccel;

		public:
			/** Initialize the arrays containing the particle data.
			*/
			virtual void init();
			
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

			void emittedParticles(FluidModel *model, const unsigned int startIndex);

			FORCE_INLINE const Real getPressure(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_pressure[fluidIndex][i];
			}

			FORCE_INLINE Real& getPressure(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_pressure[fluidIndex][i];
			}

			FORCE_INLINE void setPressure(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_pressure[fluidIndex][i] = p;
			}
			
			FORCE_INLINE Vector3r &getPressureAccel(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_pressureAccel[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r &getPressureAccel(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_pressureAccel[fluidIndex][i];
			}

			FORCE_INLINE void setPressureAccel(const unsigned int fluidIndex, const unsigned int i, const Vector3r &val)
			{
				m_pressureAccel[fluidIndex][i] = val;
			}

	};
}

#endif