#ifndef __SimulationDataDFSPH_h__
#define __SimulationDataDFSPH_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Divergence-free Smoothed Particle Hydrodynamics introduced
	* by Bender and Koschier \cite Bender:2015, \cite Bender2017.
	*/
	class SimulationDataDFSPH
	{
		public:
			SimulationDataDFSPH();
			virtual ~SimulationDataDFSPH();

		protected:	
			/** \brief factor \f$\alpha_i\f$ \cite Bender:2015 */
			std::vector<std::vector<Real>> m_factor;
			/** \brief stores \f$\kappa\f$ value of last time step for a warm start of the pressure solver */
			std::vector<std::vector<Real>> m_kappa;
			/** \brief stores \f$\kappa^v\f$ value of last time step for a warm start of the divergence solver */
			std::vector<std::vector<Real>> m_kappaV;
			/** \brief advected density */
			std::vector<std::vector<Real>> m_density_adv;

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

			FORCE_INLINE const Real getFactor(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_factor[fluidIndex][i];
			}

			FORCE_INLINE Real& getFactor(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_factor[fluidIndex][i];
			}

			FORCE_INLINE void setFactor(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_factor[fluidIndex][i] = p;
			}

			FORCE_INLINE const Real getKappa(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_kappa[fluidIndex][i];
			}

			FORCE_INLINE Real& getKappa(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_kappa[fluidIndex][i];
			}

			FORCE_INLINE void setKappa(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_kappa[fluidIndex][i] = p;
			}

			FORCE_INLINE const Real getKappaV(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_kappaV[fluidIndex][i];
			}

			FORCE_INLINE Real& getKappaV(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_kappaV[fluidIndex][i];
			}

			FORCE_INLINE void setKappaV(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_kappaV[fluidIndex][i] = p;
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
	};
}

#endif