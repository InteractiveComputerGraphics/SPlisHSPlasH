#include "VolumeSampling.h"
#include "SDFFunctions.h"
#include "Utilities/Timing.h"

using namespace Eigen;
using namespace std;
using namespace Utilities;


void VolumeSampling::sampleMesh(const unsigned int numVertices, const Vector3r *vertices,
	const unsigned int numFaces, const unsigned int *faces,
	const Real radius, const AlignedBox3r *region, 
	const std::array<unsigned int, 3> &resolution, const bool invert, 
	const unsigned int sampleMode, 
	std::vector<Vector3r> &samples)
 {
	AlignedBox3r bbox = SDFFunctions::computeBoundingBox(numVertices, vertices);

	if (region)
	{
		for (unsigned int i = 0; i < 3; i++)
		{
			bbox.min()[i] = std::max(region->min()[i], bbox.min()[i]);
			bbox.max()[i] = std::min(region->max()[i], bbox.max()[i]);
		}
	}

	Discregrid::CubicLagrangeDiscreteGrid* sdf = SDFFunctions::generateSDF(numVertices, vertices, numFaces, faces, bbox, resolution, invert);

	const Real diameter = static_cast<Real>(2.0) * radius;

	// sample object
	const unsigned int numberOfSamplePoints = (((unsigned int)((1.0f / diameter) * (bbox.max()[2] - bbox.min()[2]))) + 1) *
		(((unsigned int)((1.0f / diameter) * (bbox.max()[1] - bbox.min()[1]))) + 1) *
		(((unsigned int)((1.0f / diameter) * (bbox.max()[0] - bbox.min()[0]))) + 1);

	unsigned int currentSample = 0;
 //	Real currentPercent = 0.0;
	int counter_x = 0;
	int counter_y = 0;
	int counter_z = 0;
	Real xshift = diameter;
	Real yshift = diameter;
	Real zshift = diameter;
	Real packRadius = radius;

	if (sampleMode == 1)
	{
		// Mode 1 is a 2D hexagonal/triangular packing in the xy-plane, stacked
		// with unscaled (cubic) spacing in z. Its density ratio vs. the cubic
		// lattice (which the particle mass is based on) is 2/sqrt(3), so only
		// the in-plane (x/y) spacing must be enlarged by sqrt(2/sqrt(3)) to
		// compensate; z is left untouched.
		const Real packingFactor = sqrt(static_cast<Real>(2.0) / sqrt(static_cast<Real>(3.0)));
		const Real packDiameter = diameter * packingFactor;
		packRadius = static_cast<Real>(0.5) * packDiameter;
		xshift = packDiameter;
		yshift = sqrt(static_cast<Real>(3.0)) * packRadius;
	}
	else if (sampleMode == 2)
	{
		// Hexagonal close packing (HCP) has a higher packing fraction (~74%) than
		// the simple cubic lattice (~52%) that the particle mass is based on. To
		// still reach rest density, the lattice spacing must be enlarged by the
		// well-known factor 2^(1/6), which is the optimal spacing factor for
		// close-packed SPH samplings.
		const Real packingFactor = pow(static_cast<Real>(2.0), static_cast<Real>(1.0) / static_cast<Real>(6.0));
		const Real packDiameter = diameter * packingFactor;
		packRadius = static_cast<Real>(0.5) * packDiameter;
		xshift = packDiameter;
		yshift = sqrt(static_cast<Real>(3.0)) * packRadius;
		zshift = sqrt(static_cast<Real>(6.0)) * packDiameter / static_cast<Real>(3.0);
	}
	for (Real z = bbox.min()[2]; z <= bbox.max()[2]; z += zshift)
	{
		for (Real y = bbox.min()[1]; y <= bbox.max()[1]; y += yshift)
		{
			for (Real x = bbox.min()[0]; x <= bbox.max()[0]; x += xshift)
			{
				Vector3r particlePosition;
				if (sampleMode == 1)
				{					
					if (counter_y % 2 == 0)
						particlePosition = Vector3r(x, y + packRadius, z + packRadius);
					else
						particlePosition = Vector3r(x + packRadius, y + packRadius, z);
				}
				else if (sampleMode == 2)
				{
					particlePosition = Vector3r(x, y + packRadius, z + packRadius);

					// Hexagonal close packing (HCP):
					// shift alternating x-columns by one particle radius (within a layer and
					// between layers) and additionally shift alternating z-layers by
					// sqrt(3)/3 * radius in y-direction, so that each layer sits in the
					// pockets of the layer below. The spacing itself is already enlarged by
					// 2^(1/6) above to compensate for the higher packing fraction of HCP.
					Vector3r shift_vec(0, 0, 0);
					if (((counter_y + counter_z) % 2) != 0)
						shift_vec[0] += packRadius;
					if ((counter_z % 2) != 0)
						shift_vec[1] += sqrt(static_cast<Real>(3.0)) / static_cast<Real>(3.0) * packRadius;
					particlePosition += shift_vec;
				}
 				else
 				{
 					// Use center of voxel
 					particlePosition = Vector3r(x + radius, y + radius, z + radius);
 				}
 
				if (SDFFunctions::distance(sdf, particlePosition, 0.0) < 0.0)
					samples.push_back(particlePosition);
 				currentSample++;
 
//  				if ((Real)currentSample / (Real)numberOfSamplePoints > currentPercent)
//  				{
//  					LOG_INFO << currentPercent * 100.0 << "%";
//  					currentPercent += 0.05;
//  				}
				counter_x++;
			}
			counter_x = 0;
			counter_y++;
		}
		counter_y = 0;
		counter_z++;
	}
// 	LOG_INFO << "100%";

	delete sdf;
}
 