#ifndef __SimulationDataICSPH_h__
#define __SimulationDataICSPH_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Implicit Compressible SPH introduced
	* by Gissler et al. [GHB+20].
	*
	* References:
	* - [GHB+20] Christoph Gissler, Andreas Henne, Stefan Band, Andreas Peer and Matthias Teschner. An Implicit Compressible SPH Solver for Snow Simulation. ACM Transactions on Graphics, 39(4). URL: https://doi.org/10.1145/3386569.3392431
	*/
	class SimulationDataICSPH
	{
		public:
			SimulationDataICSPH();
			virtual ~SimulationDataICSPH();

		protected:	
			std::vector<std::vector<Real>> m_aii;
			std::vector<std::vector<Real>> m_density_adv;
			std::vector<std::vector<Real>> m_pressure;
			std::vector<std::vector<Vector3r>> m_pressureGradient;
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

			FORCE_INLINE const Real getAii(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_aii[fluidIndex][i];
			}

			FORCE_INLINE Real& getAii(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_aii[fluidIndex][i];
			}

			FORCE_INLINE void setAii(const unsigned int fluidIndex, const unsigned int i, const Real aii)
			{
				m_aii[fluidIndex][i] = aii;
			}

			FORCE_INLINE const Real getDensityAdv(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_density_adv[fluidIndex][i];
			}

			FORCE_INLINE Real& getDensityAdv(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_density_adv[fluidIndex][i];
			}

			FORCE_INLINE void setDensityAdv(const unsigned int fluidIndex, const unsigned int i, const Real d)
			{
				m_density_adv[fluidIndex][i] = d;
			}

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

			FORCE_INLINE Vector3r& getPressureGradient(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_pressureGradient[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r& getPressureGradient(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_pressureGradient[fluidIndex][i];
			}

			FORCE_INLINE void setPressureGradient(const unsigned int fluidIndex, const unsigned int i, const Vector3r& val)
			{
				m_pressureGradient[fluidIndex][i] = val;
			}

	};
}

#endif