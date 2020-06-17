#ifndef __VorticityConfinement_h__
#define __VorticityConfinement_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "VorticityBase.h"

namespace SPH
{
	/** \brief This class implements the vorticity confinement method introduced
	* by Macklin and Mueller [MM13].
	*
	* References:
	* - [MM13] Miles Macklin and Matthias Müller. Position based fluids. ACM Trans. Graph., 32(4):104:1-104:12, July 2013. URL: http://doi.acm.org/10.1145/2461912.2461984
	*/
	class VorticityConfinement : public VorticityBase
	{
	protected:
		std::vector<Vector3r> m_omega;
		std::vector<Real> m_normOmega;

	public:
		VorticityConfinement(FluidModel *model);
		virtual ~VorticityConfinement(void);

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		FORCE_INLINE const Vector3r& getAngularVelocity(const unsigned int i) const
		{
			return m_omega[i];
		}

		FORCE_INLINE Vector3r& getAngularVelocity(const unsigned int i)
		{
			return m_omega[i];
		}

		FORCE_INLINE void setAngularVelocity(const unsigned int i, const Vector3r& val)
		{
			m_omega[i] = val;
		}
	};
}

#endif
