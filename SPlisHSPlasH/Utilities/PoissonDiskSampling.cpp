#include "PoissonDiskSampling.h"

#include <algorithm>
#include <limits>

#include <iostream>
#include <fstream>
#include <string>

#define _USE_MATH_DEFINES
#include "math.h"

using namespace std;
using namespace Eigen;
using namespace SPH;

PoissonDiskSampling::PoissonDiskSampling()
{
}

void PoissonDiskSampling::sampleMesh(const unsigned int numVertices, const Vector3r *vertices, const unsigned int numFaces, const unsigned int *faces,
	const Real minRadius, const unsigned int numTrials,
	unsigned int distanceNorm, std::vector<Vector3r> &samples)
{
	m_r = minRadius;
	m_numTrials = numTrials;
	m_distanceNorm = distanceNorm;

	m_cellSize = m_r / sqrt(static_cast<Real>(3.0));

	// Init sampling
	m_maxArea = numeric_limits<Real>::min();
	determineMinX(numVertices, vertices);

	determineTriangleAreas(numVertices, vertices, numFaces, faces);

	const Real circleArea = static_cast<Real>(M_PI) * minRadius * minRadius;
	const unsigned int numInitialPoints = (unsigned int) (40.0 * (m_totalArea / circleArea)); 
	//cout << "# Initial points: " << numInitialPoints << endl;

	m_initialInfoVec.resize(numInitialPoints);
	m_phaseGroups.resize(27);

	computeFaceNormals(numVertices, vertices, numFaces, faces);

	// Generate initial set of candidate points
	generateInitialPointSet(numVertices, vertices, numFaces, faces);

	// Find minimal coordinates of object

	// Calculate CellIndices
	const Real factor = static_cast<Real>(1.0) / m_cellSize;

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < (int)m_initialInfoVec.size(); i++)
	{
		const Vector3r& v = m_initialInfoVec[i].pos;
		const int cellPos1 = PoissonDiskSampling::floor((v.x() - m_minVec[0]) * factor) + 1;
		const int cellPos2 = PoissonDiskSampling::floor((v.y() - m_minVec[1]) * factor) + 1;
		const int cellPos3 = PoissonDiskSampling::floor((v.z() - m_minVec[2]) * factor) + 1;
		m_initialInfoVec[i].cP = CellPos(cellPos1, cellPos2, cellPos3);
	}

	// Sort Initial points for CellID
	quickSort(0, (int)m_initialInfoVec.size() - 1);

	// PoissonSampling
	parallelUniformSurfaceSampling(samples);

	// release data
	m_initialInfoVec.clear();
	for (int i = 0; i < m_phaseGroups.size(); i++)
	{
		m_phaseGroups[i].clear();
	}
	m_phaseGroups.clear();
}


void PoissonDiskSampling::determineTriangleAreas(const unsigned int numVertices, const Vector3r *vertices, const unsigned int numFaces, const unsigned int *faces)
{
	m_areas.resize(numFaces);
	Real totalArea = 0.0;
	Real tmpMaxArea = numeric_limits<Real>::min();

	#pragma omp parallel default(shared)
	{
		// Compute area of each triangle
		#pragma omp for reduction(+:totalArea) schedule(static) 
		for (int i = 0; i < (int)numFaces; i++)
		{
			const Vector3r &a = vertices[faces[3 * i]];
			const Vector3r &b = vertices[faces[3 * i + 1]];
			const Vector3r &c = vertices[faces[3 * i + 2]];

			const Vector3r d1 = b - a;
			const Vector3r d2 = c - a;

			const Real area = (d1.cross(d2)).norm() / static_cast<Real>(2.0);
			m_areas[i] = area;
			totalArea += area;
			//tmpMaxArea = max(area, tmpMaxArea);

			if (area > tmpMaxArea)
			{
				#pragma omp critical
				{
					tmpMaxArea = max(area, tmpMaxArea);
				}
			}

		}
	}
	m_maxArea = max(tmpMaxArea, m_maxArea);
	m_totalArea = totalArea;
}

void PoissonDiskSampling::generateInitialPointSet(const unsigned int numVertices, const Vector3r *vertices, const unsigned int numFaces, const unsigned int *faces)
{
	m_totalArea = 0.0;
	
	// Drawing random barycentric coordinates
	std::uniform_real_distribution<Real> distribution(0.0, 1.0);
	
	random_device r;
	std::vector<std::default_random_engine> generators;

	#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	#else
	const int maxThreads = 1;
	#endif

    for (int i = 0; i < maxThreads; ++i) 
        generators.emplace_back(default_random_engine(r()));
	
	#pragma omp parallel default(shared)
	{	
		// Generating the surface points
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)m_initialInfoVec.size(); i++)
		{
			#ifdef _OPENMP
			int tid = omp_get_thread_num();
			#else
			int tid = 0;
			#endif
			// Get the generator based on thread id
			std::default_random_engine& generator = generators[tid];
			
			Real rn1 = sqrt(distribution(generator));
			Real bc1 = static_cast<Real>(1.0) - rn1;
			Real bc2 = distribution(generator)*rn1;
			Real bc3 = static_cast<Real>(1.0) - bc1 - bc2;

			// Triangle selection with probability proportional to area
			const unsigned int randIndex = getAreaIndex(m_areas, m_totalArea, generator, distribution);

			// Calculating point coordinates
			const Vector3r &v1 = vertices[faces[3 * randIndex]];
			const Vector3r &v2 = vertices[faces[3 * randIndex + 1]];
			const Vector3r &v3 = vertices[faces[3 * randIndex + 2]];

			m_initialInfoVec[i].pos = bc1*v1 + bc2*v2 + bc3*v3;
			m_initialInfoVec[i].ID = randIndex;
		}
	}
}


unsigned int PoissonDiskSampling::getAreaIndex(const vector<Real>& areas, const Real totalArea, std::default_random_engine &generator, std::uniform_real_distribution<Real> &distribution)
{
	// see https://en.wikipedia.org/wiki/Fitness_proportionate_selection
	//// Linear Version O(n)
	//Real rn = m_uniform_distribution1(m_generator)*totalArea;

	//for (unsigned int i = 0; i<areas.size(); i++)
	//{
	//	rn -= areas[i];
	//	if (rn <= 0)
	//		return i;
	//}
	//return (int)areas.size() - 1;

	// Stochastic acceptance version O(1)
	bool notaccepted = true;
	unsigned int index = 0;
	while (notaccepted)
	{
		index = (int)((Real)areas.size()*distribution(generator));
		if (distribution(generator)<areas[index] / m_maxArea)
			notaccepted = false;
	}
	return index;
}


void PoissonDiskSampling::parallelUniformSurfaceSampling(std::vector<Vector3r> &samples)
{

	// Sort initial points into HashMap storing only the index of the first point of cell
	// and build phase groups
	unordered_map<CellPos, HashEntry, CellPosHasher> hMap(2 * m_initialInfoVec.size());
	samples.clear();
	samples.reserve(m_initialInfoVec.size());

	// Already insert first Initial point as start of first cell in hashmap
	{
		const CellPos& cell = m_initialInfoVec[0].cP;
		HashEntry &entry = hMap[cell];
		entry.startIndex = 0;
		entry.samples.reserve(5);
		int index = cell[0] % 3 + 3 * (cell[1] % 3) + 9 * (cell[2] % 3);
		m_phaseGroups[index].push_back(cell);
	}

	for (int i = 1; i < (int)m_initialInfoVec.size(); i++)
	{
		const CellPos& cell = m_initialInfoVec[i].cP;
		if (cell != m_initialInfoVec[i - 1].cP)
		{
			HashEntry &entry = hMap[cell];
			entry.startIndex = i;
			entry.samples.reserve(5);
			int index = cell[0] % 3 + 3 * (cell[1] % 3) + 9 * (cell[2] % 3);
			m_phaseGroups[index].push_back(cell);
		}
	}
	// Loop over number of tries to find a sample in a cell
	for (int k = 0; k < (int)m_numTrials; k++)
	{
		// Loop over the 27 cell groups
		for (int pg = 0; pg < m_phaseGroups.size(); pg++)
		{
			const vector<CellPos>& cells = m_phaseGroups[pg];
			// Loop over the cells in each cell group
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < (int)cells.size(); i++)
			{
				const auto entryIt = hMap.find(cells[i]);
				// Check if cell exists
				if (entryIt != hMap.end())
				{
					// Check if max Index is not exceeded
					HashEntry& entry = entryIt->second;
					if (entry.startIndex + k < m_initialInfoVec.size())
					{
						if (m_initialInfoVec[entry.startIndex].cP == m_initialInfoVec[entry.startIndex + k].cP)
						{
							// choose kth point from cell
							const InitialPointInfo& test = m_initialInfoVec[entry.startIndex + k];
							// Assign sample
							if (!nbhConflict(hMap, test))
							{
								const int index = entry.startIndex + k;
								#pragma omp critical
								{
									entry.samples.push_back(index);
									samples.push_back(m_initialInfoVec[index].pos);
								}
							}
						}
					}
				}
			}
		}
	}
}

bool PoissonDiskSampling::nbhConflict(const unordered_map<CellPos, HashEntry, CellPosHasher>& hMap, const InitialPointInfo& iPI)
{
	CellPos nbPos = iPI.cP;

	// check neighboring cells inside to outside
	if (checkCell(hMap, nbPos, iPI)) 
		return true;
	for (int level = 1; level < 3; level++)
	{
		for (int ud = -level; ud < level + 1; ud += 2 * level)
		{
			for (int i = -level+1; i < level ; i++)
			{
				for (int j = -level + 1; j < level ; j++)
				{
					nbPos = CellPos(i, ud, j) + iPI.cP;
					if (checkCell(hMap, nbPos, iPI)) 
						return true;
				}
			}

			for (int i = -level; i < level + 1; i++)
			{
				for (int j = -level + 1; j < level ; j++)
				{
					nbPos = CellPos(j, i, ud) + iPI.cP;
					if (checkCell(hMap, nbPos, iPI)) 
						return true;
				}

				for (int j = -level; j < level + 1; j++)
				{
					nbPos = CellPos(ud, i, j) + iPI.cP;
					if (checkCell(hMap, nbPos, iPI)) 
						return true;
				}
			}
		}
	}
	return false;
}

bool PoissonDiskSampling::checkCell(const unordered_map<CellPos, HashEntry, CellPosHasher>& hMap, const CellPos& cell, const InitialPointInfo& iPI)
{
	const auto nbEntryIt = hMap.find(cell);
	if (nbEntryIt != hMap.end())
	{
		const HashEntry& nbEntry = nbEntryIt->second;
		for (unsigned int i = 0; i < nbEntry.samples.size(); i++)
		{
			const InitialPointInfo &info = m_initialInfoVec[nbEntry.samples[i]];
			Real dist;
			if (m_distanceNorm == 0 || iPI.ID == info.ID)
			{
				dist = (iPI.pos - info.pos).norm();
			}
			else if (m_distanceNorm == 1)
			{
				Vector3r v = (info.pos - iPI.pos).normalized();
				Real c1 = m_faceNormals[iPI.ID].dot(v);
				Real c2 = m_faceNormals[info.ID].dot(v);

				dist = (iPI.pos - info.pos).norm();
				if (fabs(c1 - c2) > 0.00001f)
					dist *= (asin(c1) - asin(c2)) / (c1 - c2);
				else
					dist /= (sqrt(static_cast<Real>(1.0) - c1*c1));
			}
			else
			{
				return true;
			}

			if (dist < m_r)
				return true;
		}
	}
	return false;
}

void PoissonDiskSampling::determineMinX(const unsigned int numVertices, const Vector3r *vertices)
{
	m_minVec = Vector3r(numeric_limits<Real>::max(), numeric_limits<Real>::max(), numeric_limits<Real>::max());

	for (int i = 0; i < (int)numVertices; i++)
	{
		const Vector3r& v = vertices[i];
		m_minVec[0] = std::min(m_minVec[0], v[0]);
		m_minVec[1] = std::min(m_minVec[1], v[1]);
		m_minVec[2] = std::min(m_minVec[2], v[2]);
	}
}

void PoissonDiskSampling::quickSort(int left, int right)
{
	if (left < right)
	{
		int index = partition(left, right);
		quickSort(left, index - 1);
		quickSort(index, right);
	}
}

int PoissonDiskSampling::partition(int left, int right)
{
	int i = left;
	int j = right;
	Vector3r tmpPos;
	CellPos tmpCell;
	InitialPointInfo tmpInfo;
	CellPos pivot = m_initialInfoVec[left + (right - left) / 2].cP;
	
	while (i <= j)
	{
		while (compareCellID(m_initialInfoVec[i].cP, pivot))
			i++;

		while (compareCellID(pivot, m_initialInfoVec[j].cP))
			j--;

		if (i <= j)
		{
			tmpInfo = m_initialInfoVec[i];
			m_initialInfoVec[i] = m_initialInfoVec[j];
			m_initialInfoVec[j] = tmpInfo;
			i++;
			j--;
		}
	}
	return i;
}

bool PoissonDiskSampling::compareCellID(CellPos& a, CellPos& b)
{
	for (unsigned int i = 0; i < 3; i++)
	{
		if (a[i] < b[i]) return true;
		if (a[i] > b[i]) return false;
	}

	return false;
}

void PoissonDiskSampling::computeFaceNormals(const unsigned int numVertices, const Vector3r *vertices, const unsigned int numFaces, const unsigned int *faces)
{
	m_faceNormals.resize(numFaces);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numFaces; i++)
		{
			// Get first three points of face
			const Vector3r &a = vertices[faces[3 * i]];
			const Vector3r &b = vertices[faces[3 * i + 1]];
			const Vector3r &c = vertices[faces[3 * i + 2]];

			// Create normal
			Vector3r v1 = b - a;
			Vector3r v2 = c - a;

			m_faceNormals[i] = v1.cross(v2);
			m_faceNormals[i].normalize();
		}
	}
}
