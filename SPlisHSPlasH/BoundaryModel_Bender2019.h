#ifndef __BoundaryModel_Bender2019_h__
#define __BoundaryModel_Bender2019_h__

#include "Common.h"
#include <vector>

#include "SPHKernels.h"
#include "Discregrid/All"
#include "BoundaryModel.h"


namespace SPH 
{	
	class TimeStep;

	/** \brief The boundary model stores the information required for boundary handling
	* using the approach of Bender et al. 2019 [BKWK19].
	*
	* References:
	* - [BKWK19] Jan Bender, Tassilo Kugelstadt, Marcel Weiler, and Dan Koschier. Volume maps: an implicit boundary representation for SPH. In Proceedings of ACM SIGGRAPH Conference on Motion, Interaction and Games, MIG '19. ACM, 2019. URL: https://dl.acm.org/doi/10.1145/3359566.3360077
	*/
	class BoundaryModel_Bender2019 : public BoundaryModel
	{
		public:
			BoundaryModel_Bender2019();
			virtual ~BoundaryModel_Bender2019();


	protected:
			// Density or volume map 
			Discregrid::DiscreteGrid *m_map;
			// values required for volume maps
			std::vector<std::vector<Real>> m_boundaryVolume;
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

			FORCE_INLINE const Real& getBoundaryVolume(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_boundaryVolume[fluidIndex][i];
			}

			FORCE_INLINE Real& getBoundaryVolume(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_boundaryVolume[fluidIndex][i];
			}

			FORCE_INLINE void setBoundaryVolume(const unsigned int fluidIndex, const unsigned int i, const Real &val)
			{
				m_boundaryVolume[fluidIndex][i] = val;
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