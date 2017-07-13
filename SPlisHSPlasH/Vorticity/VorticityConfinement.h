#ifndef __VorticityConfinement_h__
#define __VorticityConfinement_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief This class implements the vorticity confinement method introduced
	* by Macklin and Mueller \cite Macklin:2013:PBF.
	*/
	class VorticityConfinement : public NonPressureForceBase
	{
	protected:
		std::vector<Vector3r> m_omega;
		std::vector<Real> m_normOmega;
		Real m_viscosityT;

	public:
		VorticityConfinement(FluidModel *model);
		virtual ~VorticityConfinement(void);

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		Real getViscosityT() const { return m_viscosityT; }
		void setViscosityT(Real val) { m_viscosityT = val; }

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
