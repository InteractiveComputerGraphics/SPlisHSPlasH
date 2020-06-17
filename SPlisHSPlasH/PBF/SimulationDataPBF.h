#ifndef __SimulationDataPBF_h__
#define __SimulationDataPBF_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Position-Based Fluids introduced
	* by Macklin and Mueller [MM13,BMO+14,BMM15].
	*
	* References:
	* - [MM13] Miles Macklin and Matthias Müller. Position based fluids. ACM Trans. Graph., 32(4):104:1-104:12, July 2013. URL: http://doi.acm.org/10.1145/2461912.2461984
	* - [BMO+14] Jan Bender, Matthias Müller, Miguel A. Otaduy, Matthias Teschner, and Miles Macklin. A survey on position-based simulation methods in computer graphics. Computer Graphics Forum, 33(6):228-251, 2014. URL: http://dx.doi.org/10.1111/cgf.12346
	* - [BMM15] Jan Bender, Matthias Müller, and Miles Macklin. Position-based simulation methods in computer graphics. In EUROGRAPHICS 2015 Tutorials. Eurographics Association, 2015. URL: http://dx.doi.org/10.2312/egt.20151045
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