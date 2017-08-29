#ifndef __SurfaceTension_Becker2007_h__
#define __SurfaceTension_Becker2007_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"

namespace SPH
{
	/** \brief This class implements the surface tension method introduced
	* by Becker and Teschner \cite Becker:2007.
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
