#ifndef __SurfaceTension_Becker2007_h__
#define __SurfaceTension_Becker2007_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief This class implements the surface tension method introduced
	* by Becker and Teschner [BT07].
	*
	* References: 
	* - [BT07] Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In ACM SIGGRAPH/Eurographics Symposium on Computer Animation, SCA '07, 209-217. Aire-la-Ville, Switzerland, Switzerland, 2007. Eurographics Association. URL: http://dl.acm.org/citation.cfm?id=1272690.1272719
	*/
	class SurfaceTension_Becker2007 : public NonPressureForceBase
	{
	protected:
		Real m_surfaceTension;
		Real m_surfaceTensionBoundary;

		virtual void initParameters();

	public:
		static std::string METHOD_NAME;
		static int SURFACE_TENSION;
		static int SURFACE_TENSION_BOUNDARY;

		SurfaceTension_Becker2007(FluidModel *model);
		virtual ~SurfaceTension_Becker2007(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_Becker2007(model); }
		virtual std::string getMethodName() { return METHOD_NAME; }

		virtual void step();
		virtual void reset();
	};
}

#endif
