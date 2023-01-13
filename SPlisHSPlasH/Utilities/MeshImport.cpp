#include "MeshImport.h"

#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include <Utilities/PLYLoader.h>
#include <Utilities/OBJLoader.h>

using namespace SPH;
using namespace Utilities;
using namespace std;


bool MeshImport::importMesh(const std::string& filename, TriangleMesh& mesh,
	const Vector3r& translation, const Matrix3r& rotation, const Vector3r& scale)
{
	if (!FileSystem::fileExists(filename))
	{
		LOG_ERR << "File not found: " << filename;
		return false;
	}
	string ext = FileSystem::getFileExt(filename);
	transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

	if (ext == "PLY")
		return importMesh_PLY(filename, mesh, translation, rotation, scale);
	else if (ext == "OBJ")
		return importMesh_OBJ(filename, mesh, translation, rotation, scale);
	else
	{
		LOG_ERR << "File " << filename << " has an unknown file type.";
		return false;
	}
}

bool MeshImport::importMesh_PLY(const std::string& filename, TriangleMesh& mesh,
	const Vector3r& translation, const Matrix3r& rotation, const Vector3r& scale)
{
	std::vector<std::array<float, 3>> x;
	std::vector<std::array<int, 3>> faces;
	std::array<float, 3> s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	PLYLoader::loadPly(filename, x, faces, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	mesh.initMesh(nPoints, nFaces);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		mesh.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		int posIndices[3];
		for (int j = 0; j < 3; j++)
			posIndices[j] = faces[i][j];

		mesh.addFace(&posIndices[0]);
	}

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;
	return true;
}


bool MeshImport::importMesh_OBJ(const std::string& filename, TriangleMesh& mesh,
	const Vector3r& translation, const Matrix3r& rotation, const Vector3r& scale)
{
	std::vector<std::array<float, 3>> x;
	std::vector<OBJLoader::Vec3f> normals;
	std::vector<MeshFaceIndices> faces;
	OBJLoader::Vec3f s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	OBJLoader::loadObj(filename, &x, &faces, &normals, nullptr, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	mesh.initMesh(nPoints, nFaces);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		mesh.addVertex(rotation * Vector3r(x[i][0], x[i][1], x[i][2]) + translation);
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		int posIndices[3];
		for (int j = 0; j < 3; j++)
		{
			posIndices[j] = faces[i].posIndices[j];
		}

		mesh.addFace(&posIndices[0]);
	}

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;

	return true;
}