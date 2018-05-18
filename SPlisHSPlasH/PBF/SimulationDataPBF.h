#ifndef __SimulationDataPBF_h__
#define __SimulationDataPBF_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Position-Based Fluids introduced
	* by Macklin and Mueller \cite Macklin:2013:PBF, \cite BMOTM2014, \cite BMM2015.
	*/
	class SimulationDataPBF
	{
		public:
			SimulationDataPBF();
			virtual ~SimulationDataPBF();

		protected:	
			std::vector<std::vector<Real>> m_lambda;		
			std::vector<std::vector<Vector3r>> m_deltaX;
			std::vector<std::vector<Vector3r>> m_oldX;
			std::vector<std::vector<Vector3r>> m_lastX;

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

			FORCE_INLINE const Real& getLambda(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_lambda[fluidIndex][i];
			}

			FORCE_INLINE Real& getLambda(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_lambda[fluidIndex][i];
			}

			FORCE_INLINE void setLambda(const unsigned int fluidIndex, const unsigned int i, const Real &val)
			{
				m_lambda[fluidIndex][i] = val;
			}

			FORCE_INLINE Vector3r &getDeltaX(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_deltaX[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r &getDeltaX(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_deltaX[fluidIndex][i];
			}

			FORCE_INLINE void setDeltaX(const unsigned int fluidIndex, const unsigned int i, const Vector3r &val)
			{
				m_deltaX[fluidIndex][i] = val;
			}

			FORCE_INLINE Vector3r &getLastPosition(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_lastX[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r &getLastPosition(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_lastX[fluidIndex][i];
			}

			FORCE_INLINE void setLastPosition(const unsigned int fluidIndex, const unsigned int i, const Vector3r &pos)
			{
				m_lastX[fluidIndex][i] = pos;
			}

			FORCE_INLINE Vector3r &getOldPosition(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_oldX[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r &getOldPosition(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_oldX[fluidIndex][i];
			}

			FORCE_INLINE void setOldPosition(const unsigned int fluidIndex, const unsigned int i, const Vector3r &pos)
			{
				m_oldX[fluidIndex][i] = pos;
			}

	};
}

#endif