#ifndef __MeshImport_h__
#define __MeshImport_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TriangleMesh.h"

namespace SPH
{
	class MeshImport
	{
	protected:
		/** Load a mesh from an OBJ file in the TriangleMesh. */
		static bool importMesh_OBJ(const std::string& filename, TriangleMesh& mesh,
			const Vector3r& translation, const Matrix3r& rotation, const Vector3r& scale);

		/** Load a mesh from an PLY file in the TriangleMesh. */
		static bool importMesh_PLY(const std::string& filename, TriangleMesh& mesh,
			const Vector3r& translation, const Matrix3r& rotation, const Vector3r& scale);

	public:
		/** Load a mesh from a file in the TriangleMesh. */
		static bool importMesh(const std::string& filename, TriangleMesh& mesh,
			const Vector3r& translation, const Matrix3r& rotation, const Vector3r& scale);
	};
}

#endif

