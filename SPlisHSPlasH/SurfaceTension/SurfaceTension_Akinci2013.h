#ifndef __SurfaceTension_Akinci2013_h__
#define __SurfaceTension_Akinci2013_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"

namespace SPH
{
	/** \brief This class implements the surface tension method introduced
	* by Akinci et al. \cite Akinci:2013.
	*/
	class SurfaceTension_Akinci2013 : public SurfaceTensionBase
	{
	protected: 
		std::vector<Vector3r> m_normals;

	public:
		SurfaceTension_Akinci2013(FluidModel *model);
		virtual ~SurfaceTension_Akinci2013(void);

		virtual void step();
		virtual void reset();

		void computeNormals();

		virtual void performNeighborhoodSearchSort();

		FORCE_INLINE Vector3r &getNormal(const unsigned int i)
		{
			return m_normals[i];
		}

		FORCE_INLINE const Vector3r &getNormal(const unsigned int i) const
		{
			return m_normals[i];
		}

		FORCE_INLINE void setNormal(const unsigned int i, const Vector3r &val)
		{
			m_normals[i] = val;
		}
	};
}

#endif
