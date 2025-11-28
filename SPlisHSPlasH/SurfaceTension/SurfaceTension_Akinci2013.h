#ifndef __SurfaceTension_Akinci2013_h__
#define __SurfaceTension_Akinci2013_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief This class implements the surface tension method introduced
	* by Akinci et al. [ATT13].
	*
	* References:
	* - [AAT13] Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for sph fluids. ACM Trans. Graph., 32(6):182:1-182:8, November 2013. URL: http://doi.acm.org/10.1145/2508363.2508395
	*/
	class SurfaceTension_Akinci2013 : public NonPressureForceBase
	{
	protected: 
		Real m_surfaceTension;
		Real m_surfaceTensionBoundary;

		std::vector<Vector3r> m_normals;

		virtual void initParameters();

	public:
		static std::string METHOD_NAME;
		static int SURFACE_TENSION;
		static int SURFACE_TENSION_BOUNDARY;

		SurfaceTension_Akinci2013(FluidModel *model);
		virtual ~SurfaceTension_Akinci2013(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_Akinci2013(model); }
		virtual std::string getMethodName() { return METHOD_NAME; }

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
