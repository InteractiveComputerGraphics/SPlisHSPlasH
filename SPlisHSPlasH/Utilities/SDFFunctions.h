#ifndef SDFFunctions_H
#define SDFFunctions_H

#include "../Common.h"
#include "Discregrid/All"

namespace Utilities
{
	/** \brief Functions for generating and querying an SDF. 
	*/
	class SDFFunctions
	{
	public:
		/** Generate SDF from mesh.
		*/
		static Discregrid::CubicLagrangeDiscreteGrid* generateSDF(const unsigned int numVertices, 
			const Vector3r *vertices, const unsigned int numFaces, const unsigned int *faces, 
			const AlignedBox3r &bbox, const std::array<unsigned int, 3> &resolution, const bool invert=false);
	
		/** Compute the bounding box of a mesh. 
		 */
		static AlignedBox3r computeBoundingBox(const unsigned int numVertices, const Vector3r *vertices);

		/** Determine distance of a point x to the surface represented by the SDF and corresponding surface normal and 
		* next point on the surface.
		*/
		static double distance(Discregrid::CubicLagrangeDiscreteGrid* sdf, const Vector3r &x,
			const Real thickness, Vector3r &normal, Vector3r &nextSurfacePoint);

		/** Determine distance of a point x to the surface represented by the SDF. 
		 */
		static double distance(Discregrid::CubicLagrangeDiscreteGrid* sdf,
			const Vector3r &x, const Real thickness);

	};
}

#endif // SDFFunctions_H
