#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TriangleMesh.h"

namespace SPH
{
	class WindingNumbers
	{
	public:
		/** Determine the winding number for a point p and a triangle abc. 
		 */
		static Real computeGeneralizedWindingNumber(const Vector3r& p, const Vector3r& a, const Vector3r& b, const Vector3r& c);

		/** Determine the winding number of a point p in a triangle mesh.
		*/
		static Real computeGeneralizedWindingNumber(const Vector3r& p, const TriangleMesh &mesh);
	};
}