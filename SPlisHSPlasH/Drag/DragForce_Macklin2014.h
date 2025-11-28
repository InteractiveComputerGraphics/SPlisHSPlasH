#ifndef __DragForce_Macklin2014_h__
#define __DragForce_Macklin2014_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief This class implements the drag force computation introduced
	* by Macklin et al. [MMCK14].
	*
	* References:
	* - [MMCK14] Miles Macklin, Matthias Müller, Nuttapong Chentanez, and Tae-Yong Kim. Unified Particle Physics for Real-Time Applications. ACM Trans. Graph., 33(4):1-12, 2014. URL: http://doi.acm.org/10.1145/2601097.2601152
	*/
	class DragForce_Macklin2014 : public NonPressureForceBase
	{
	protected: 
		Real m_dragCoefficient;

		virtual void initParameters();

	public:
		static std::string METHOD_NAME;
		static int DRAG_COEFFICIENT;

		DragForce_Macklin2014(FluidModel *model);
		virtual ~DragForce_Macklin2014(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new DragForce_Macklin2014(model); }
		virtual std::string getMethodName() { return METHOD_NAME; }

		virtual void step();
		virtual void reset();
	};
}

#endif
