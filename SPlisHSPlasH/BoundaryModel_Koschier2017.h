#ifndef __BoundaryModel_Koschier2017_h__
#define __BoundaryModel_Koschier2017_h__

#include "Common.h"
#include <vector>

#include "SPHKernels.h"
#include "Discregrid/All"
#include "BoundaryModel.h"


namespace SPH 
{	
	class TimeStep;

	/** \brief The boundary model stores the information required for boundary handling
	* using the approach of Koschier and Bender 2017 \cite KB17.
	*/
	class BoundaryModel_Koschier2017 : public BoundaryModel
	{
		public:
			BoundaryModel_Koschier2017();
			virtual ~BoundaryModel_Koschier2017();

		protected:
			// Density map 
			Discregrid::DiscreteGrid *m_map;
			// values required for density maps
			std::vector<std::vector<Real>> m_boundaryDensity;
			std::vector<std::vector<Vector3r>> m_boundaryDensityGradient;
			std::vector<std::vector<Vector3r>> m_boundaryXj;
			// maxmimal distance of a mesh point to the center of mass (required for CFL)
			Real m_maxDist;
			Real m_maxVel;

		public:
			void initModel(RigidBodyObject *rbo);

			virtual void reset();

			Discregrid::DiscreteGrid *getMap() { return m_map; }
			void setMap(Discregrid::DiscreteGrid *map) { m_map = map; }

			Real getMaxDist() const { return m_maxDist; }
			void setMaxDist(Real val) { m_maxDist = val; }

			Real getMaxVel() const { return m_maxVel; }
			void setMaxVel(Real val) { m_maxVel = val; }

			FORCE_INLINE const Real& getBoundaryDensity(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_boundaryDensity[fluidIndex][i];
			}

			FORCE_INLINE Real& getBoundaryDensity(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_boundaryDensity[fluidIndex][i];
			}

			FORCE_INLINE void setBoundaryDensity(const unsigned int fluidIndex, const unsigned int i, const Real &val)
			{
				m_boundaryDensity[fluidIndex][i] = val;
			}

			FORCE_INLINE Vector3r &getBoundaryDensityGradient(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_boundaryDensityGradient[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r &getBoundaryDensityGradient(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_boundaryDensityGradient[fluidIndex][i];
			}

			FORCE_INLINE void setBoundaryDensityGradient(const unsigned int fluidIndex, const unsigned int i, const Vector3r &val)
			{
				m_boundaryDensityGradient[fluidIndex][i] = val;
			}

			FORCE_INLINE Vector3r &getBoundaryXj(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_boundaryXj[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r &getBoundaryXj(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_boundaryXj[fluidIndex][i];
			}

			FORCE_INLINE void setBoundaryXj(const unsigned int fluidIndex, const unsigned int i, const Vector3r &val)
			{
				m_boundaryXj[fluidIndex][i] = val;
			}
	};
}

#endif