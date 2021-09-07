#include "SPlisHSPlasH/Common.h"
#include "SPHSamplingBase.h"

#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/Version.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "Utilities/FileSystem.h"

#define USE_ADHESION


using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;

SPHSamplingBase::SPHSamplingBase() :
	SamplingBase()
{
	m_viscosity = 0.25;
	m_cohesion = 4.0;
	m_adhesion = 2.0;
	m_cflFactor = 0.5;
	m_steps = 200;
	m_neighborhoodSearch = nullptr;
	m_counter = 0;
}

SPHSamplingBase::~SPHSamplingBase()
{

}

void SPHSamplingBase::writeParticleDataVTK(const std::string& fileName)
{
	const unsigned int numParticles = (int)m_x.size();
	if (0 == numParticles)
		return;

#ifdef USE_DOUBLE
	const char* real_str = " double\n";
#else 
	const char* real_str = " float\n";
#endif

	// Open the file
	std::ofstream outfile{ fileName, std::ios::binary };
	if (!outfile.is_open())
	{
		LOG_WARN << "Cannot open a file to save VTK particles.";
		return;
	}

	outfile << "# vtk DataFile Version 4.1\n";
	outfile << "SPlisHSPlasH particle data\n"; // title of the data set, (any string up to 256 characters+\n)
	outfile << "BINARY\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";

	//////////////////////////////////////////////////////////////////////////
	// export position attribute as POINTS
	{
		std::vector<Vector3r> positions;
		positions.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			positions.emplace_back(m_x[i]);
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			for (unsigned int c = 0; c < 3; c++)
				swapByteOrder(&positions[i][c]);
		// export to vtk
		outfile << "POINTS " << numParticles << real_str;
		outfile.write(reinterpret_cast<char*>(positions[0].data()), 3 * numParticles * sizeof(Real));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// export particle IDs as CELLS
	{
		std::vector<Eigen::Vector2i> cells;
		cells.reserve(numParticles);
		unsigned int nodes_per_cell_swapped = 1;
		swapByteOrder(&nodes_per_cell_swapped);
		for (unsigned int i = 0u; i < numParticles; i++)
		{
			unsigned int idSwapped = i;
			swapByteOrder(&idSwapped);
			cells.emplace_back(nodes_per_cell_swapped, idSwapped);
		}

		// particles are cells with one element and the index of the particle
		outfile << "CELLS " << numParticles << " " << 2 * numParticles << "\n";
		outfile.write(reinterpret_cast<char*>(cells[0].data()), 2 * numParticles * sizeof(unsigned int));
		outfile << "\n";
	}
	//////////////////////////////////////////////////////////////////////////
	// export cell types
	{
		// the type of a particle cell is always 1
		std::vector<int> cellTypes;
		int cellTypeSwapped = 1;
		swapByteOrder(&cellTypeSwapped);
		cellTypes.resize(numParticles, cellTypeSwapped);
		outfile << "CELL_TYPES " << numParticles << "\n";
		outfile.write(reinterpret_cast<char*>(cellTypes.data()), numParticles * sizeof(int));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// write additional attributes as per-particle data
	{
		outfile << "POINT_DATA " << numParticles << "\n";
		// write IDs
		outfile << "SCALARS id unsigned_int 1\n";
		outfile << "LOOKUP_TABLE id_table\n";
		// copy data
		std::vector<unsigned int> attrData;
		attrData.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			attrData.emplace_back(i);
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			swapByteOrder(&attrData[i]);
		// export to vtk
		outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(unsigned int));
		outfile << "\n";
	}

	//////////////////////////////////////////////////////////////////////////
	// per point fields (all attributes except for positions)
	const auto numFields = 2;
	outfile << "FIELD FieldData " << std::to_string(numFields) << "\n";

	// Density
	{
		// write header information
		outfile << "density" << " 1 " << numParticles << real_str;

		// copy data
		std::vector<Real> attrData;
		attrData.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			attrData.emplace_back(m_densities[i]);
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			swapByteOrder(&attrData[i]);
		// export to vtk
		outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(Real));
	}

	// Corrections
	{
		// write header information
		outfile << "deltaX" << " 3 " << numParticles << real_str;

		// copy data
		std::vector<Vector3r> attrData;
		attrData.reserve(numParticles);
		for (unsigned int i = 0u; i < numParticles; i++)
			attrData.emplace_back(m_corrs[i]);
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			swapByteOrder(&attrData[i]);
		// export to vtk
		outfile.write(reinterpret_cast<char*>(attrData[0].data()), 3 * numParticles * sizeof(Real));
	}

	outfile.close();
}

void SPHSamplingBase::computeDensities(std::vector<Real> &densities, const Real mass)
{
 	const Real density0 = 1000.0;
 	const int numParticles = (int) m_x.size();
	#pragma omp parallel default(shared)
 	{
 		#pragma omp for schedule(static)  
 		for (int i = 0; i < numParticles; i++)
 		{
 			Real &density = densities[i];
 
			// Compute current density for particle i
 			density = m_W_zero;
 			const Vector3r &xi = m_x[i];
 
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r &xj = m_x[neighborIndex];
				density += m_kernelFct(xi - xj);
			}

			density *= mass;
 		}
 	}
}

void SPHSamplingBase::computeDFSPHFactor()
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	const int numParticles = (int) m_x.size();

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure stiffness denominator
		//////////////////////////////////////////////////////////////////////////

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			//////////////////////////////////////////////////////////////////////////
			// Compute gradient dp_i/dx_j * (1/k)  and dp_j/dx_j * (1/k)
			//////////////////////////////////////////////////////////////////////////
			const Vector3r &xi = m_x[i];
			Real sum_grad_p_k = 0.0;
			Vector3r grad_p_i;
			grad_p_i.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r& xj = m_x[neighborIndex];
				const Vector3r grad_p_j = -m_volume * m_gradKernelFct(xi - xj);
				sum_grad_p_k += grad_p_j.squaredNorm();
				grad_p_i -= grad_p_j;
			}

			sum_grad_p_k += grad_p_i.squaredNorm();

			//////////////////////////////////////////////////////////////////////////
			// Compute pressure stiffness denominator
			//////////////////////////////////////////////////////////////////////////
			Real &factor = m_factors[i];
			if (sum_grad_p_k > m_eps)
				factor = -static_cast<Real>(1.0) / (sum_grad_p_k);
			else
				factor = 0.0;
		}
	}
}


void SPHSamplingBase::pressureSolve()
{
	unsigned int iterations = 0;
	unsigned int maxIterations = 100;
	Real maxError = 0.01;

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////
	
	Real avg_density_err = 0.0;
	bool chk = false;

	
	while ((!chk || (iterations < 1)) && (iterations < maxIterations))
	{
		chk = true;
		avg_density_err = 0.0;

		pressureSolveIteration(avg_density_err);

		// Maximal allowed density fluctuation
		const Real eta = maxError * static_cast<Real>(0.01) * m_density0;  // maxError is given in percent
		chk = chk && (avg_density_err <= eta);

		iterations++;
	}
}


void SPHSamplingBase::pressureSolveIteration(Real &avg_density_err)
{
	const int numParticles = (int) m_x.size();
	if (numParticles == 0)
		return;

	Real density_error = 0.0;

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure forces
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			//////////////////////////////////////////////////////////////////////////
			// Evaluate rhs
			//////////////////////////////////////////////////////////////////////////
			const Real b_i = m_densities[i] / m_density0 - static_cast<Real>(1.0);
			const Real ki = b_i * m_factors[i];

			Vector3r& corr_i = m_corrs[i];
			const Vector3r& xi = m_x[i];

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r& xj = m_x[neighborIndex];
				const Vector3r grad_p_j = -m_volume * m_gradKernelFct(xi - xj);

				//const Real b_j = std::max(densities[neighborIndex] / density0 - static_cast<Real>(1.0), static_cast<Real>(0.0)); 
				const Real b_j = m_densities[neighborIndex] / m_density0 - static_cast<Real>(1.0);
				const Real kj = b_j * m_factors[neighborIndex];
				const Real kSum = ki + kj;
				if (fabs(kSum) > m_eps)
				{
					const Vector3r grad_p_j = -m_volume * m_gradKernelFct(xi - xj);

					corr_i -= kSum * grad_p_j;			// ki, kj already contain inverse density						
				}
			}
		}
	}
}
