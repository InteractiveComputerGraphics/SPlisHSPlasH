#ifndef __SurfaceTension_Becker2007_h__
#define __SurfaceTension_Becker2007_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"

namespace SPH
{
	/** \brief This class implements the surface tension method introduced
	* by Becker and Teschner [BT07].
	*
	* References: 
	* - [BT07] Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In ACM SIGGRAPH/Eurographics Symposium on Computer Animation, SCA '07, 209-217. Aire-la-Ville, Switzerland, Switzerland, 2007. Eurographics Association. URL: http://dl.acm.org/citation.cfm?id=1272690.1272719
	*/
	class SurfaceTension_Becker2007 : public SurfaceTensionBase
	{
	public:
		SurfaceTension_Becker2007(FluidModel *model);
		virtual ~SurfaceTension_Becker2007(void);

		virtual void step();
		virtual void reset();
	};
}

#endif
