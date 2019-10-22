#ifndef SURFACESAMPLING_H
#define SURFACESAMPLING_H

#include "SPlisHSPlasH/Utilities/PoissonDiskSampling.h"
#include "SPlisHSPlasH/Utilities/RegularTriangleSampling.h"
#include "SPlisHSPlasH/Utilities/RegularSampling2D.h"

namespace SPH
{
	enum SurfaceSamplingMode { PoissonDisk, RegularTriangle, Regular2D };
}

#endif // SURFACESAMPLING_H
