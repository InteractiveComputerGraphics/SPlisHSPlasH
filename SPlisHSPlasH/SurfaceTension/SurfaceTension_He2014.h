#ifndef __SurfaceTension_He2014_h__
#define __SurfaceTension_He2014_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"

namespace SPH
{
	/** \brief This class implements the surface tension method introduced
	* by He et al. \cite He:2014.
	*/
	class SurfaceTension_He2014 : public SurfaceTensionBase
	{
	protected: 
		std::vector<Real> m_color;
		std::vector<Real> m_gradC2;

	public:
		SurfaceTension_He2014(FluidModel *model);
		virtual ~SurfaceTension_He2014(void);

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		FORCE_INLINE const Real getColor(const unsigned int i) const
		{
			return m_color[i];
		}

		FORCE_INLINE Real& getColor(const unsigned int i)
		{
			return m_color[i];
		}

		FORCE_INLINE void setColor(const unsigned int i, const Real p)
		{
			m_color[i] = p;
		}

		FORCE_INLINE const Real getGradC2(const unsigned int i) const
		{
			return m_gradC2[i];
		}

		FORCE_INLINE Real& getGradC2(const unsigned int i)
		{
			return m_gradC2[i];
		}

		FORCE_INLINE void setGradC2(const unsigned int i, const Real p)
		{
			m_gradC2[i] = p;
		}
	};
}

#endif
