#ifndef __VorticityRefinement_Liu2021_h__
#define __VorticityRefinement_Liu2021_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief This class implements the vorticity refinement model introduced
	* by Liu et al. [LWB*21].
	*
	* References:
	* - [LWB*21] Liu, S., Wang, X., Ban, X., Xu, Y., Zhou, J., Kosinka, J. and Telea, A.C. (2021), Turbulent Details Simulation for SPH Fluids via Vorticity Refinement. Computer Graphics Forum, 40: 54-67. https://doi.org/10.1111/cgf.14095
	*/
	class VorticityRefinement_Liu2021 : public NonPressureForceBase
	{
	protected:
		std::vector<Vector3r> m_vorticity;			// vorticity at end of last time step
		std::vector<Vector3r> m_dissipatedVorticity;
		std::vector<Vector3r> m_streamFunction;

		Real m_vorticityCoeff;
		Real m_viscosityCoefficient;

		virtual void initParameters();

	public:
		static std::string METHOD_NAME;
		static int VORTICITY_COEFFICIENT;
		static int KINEMATIC_VISCOSITY;

		VorticityRefinement_Liu2021(FluidModel *model);
		virtual ~VorticityRefinement_Liu2021(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new VorticityRefinement_Liu2021(model); }
		virtual std::string getMethodName() { return METHOD_NAME; }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		FORCE_INLINE const Vector3r& getVorticity(const unsigned int i) const
		{
			return m_vorticity[i];
		}

		FORCE_INLINE Vector3r& getVorticity(const unsigned int i)
		{
			return m_vorticity[i];
		}

		FORCE_INLINE void setVorticity(const unsigned int i, const Vector3r& val)
		{
			m_vorticity[i] = val;
		}

		FORCE_INLINE const Vector3r& getDissipatedVorticity(const unsigned int i) const
		{
			return m_dissipatedVorticity[i];
		}

		FORCE_INLINE Vector3r& getDissipatedVorticity(const unsigned int i)
		{
			return m_dissipatedVorticity[i];
		}

		FORCE_INLINE void setDissipatedVorticity(const unsigned int i, const Vector3r& val)
		{
			m_dissipatedVorticity[i] = val;
		}

		FORCE_INLINE const Vector3r& getStreamFunction(const unsigned int i) const
		{
			return m_streamFunction[i];
		}

		FORCE_INLINE Vector3r& getStreamFunction(const unsigned int i)
		{
			return m_streamFunction[i];
		}

		FORCE_INLINE void setStreamFunction(const unsigned int i, const Vector3r& val)
		{
			m_streamFunction[i] = val;
		}
	};
}

#endif
