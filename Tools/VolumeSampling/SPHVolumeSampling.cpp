#include "SPlisHSPlasH/Common.h"
#include "SPHVolumeSampling.h"

#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/Version.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "Utilities/FileSystem.h"


using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;

SPHVolumeSampling::SPHVolumeSampling() :
	SPHSamplingBase()
{
}

SPHVolumeSampling::~SPHVolumeSampling()
{

}

void SPHVolumeSampling::generateSamples()
{
	m_volume = static_cast<Real>(1.0) * m_diameter * m_diameter * m_diameter;
	m_mass = m_volume * m_density0;

	const Real oldRadius = m_radius;

	// use smaller radius for sampling to generate more particles		
	m_radius = static_cast<Real>(0.95) * m_radius;
	m_diameter = static_cast<Real>(2.0) * m_radius;
	sampleObject(1);
	m_radius = oldRadius;
	m_diameter = static_cast<Real>(2.0) * m_radius;

	LOG_INFO << "Generated samples: " << m_x.size();

	initSPHOptimization();

	std::string exePath = FileSystem::getProgramPath();
	std::string exportPath = FileSystem::normalizePath(m_outputPath + "/SPHVolumeSampling_Kugelstadt2021");
	FileSystem::makeDirs(exportPath);

	Real maxError = 0.001;

	std::string fileName = "SamplingStep0";
	std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName + ".bgeo");
	PartioReaderWriter::writeParticles(exportFileName, (unsigned int)m_x.size(), m_x.data(), NULL, 0.0);

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////

	Real avg_density_err = 0.0;
	for (m_counter = 0; m_counter < m_steps; m_counter++)
	{
		avg_density_err = 0.0;
		START_TIMING("timeStep")
		step(avg_density_err);
		STOP_TIMING_AVG;
		std::cout << "\r" << "SPH simulation steps " << std::setw(23) << m_counter + 1;

		std::string fileName = "SamplingStep";
		fileName = fileName + std::to_string(m_counter + 1);
		if (m_output_format == 0)
		{
			std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName + ".bgeo");
			PartioReaderWriter::writeParticles(exportFileName, (unsigned int)m_x.size(), m_x.data(), NULL, 0.0);
		}
		else
		{
			std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName + ".vtk");
			writeParticleDataVTK(exportFileName);
		}

		// Maximal allowed density fluctuation
		const Real eta = maxError * static_cast<Real>(0.01) * m_density0;  // maxError is given in percent
	}

	std::cout << std::endl;
	delete m_neighborhoodSearch;
}

void SPHVolumeSampling::initSPHOptimization()
{
	m_mass = static_cast<Real>(1.0)* m_diameter* m_diameter* m_diameter* m_density0;

	// Init neighborhood search
	const Real supportRadius = static_cast<Real>(4.0) * m_radius;
	m_neighborhoodSearch = new NeighborhoodSearch(supportRadius, false);
	m_neighborhoodSearch->set_radius(supportRadius);
	m_neighborhoodSearch->add_point_set(&m_x[0][0], m_x.size(), true, true);

	CubicKernel::setRadius(supportRadius);
	CohesionKernel::setRadius(supportRadius);

	m_W_zero = CubicKernel::W_zero();
	m_kernelFct = CubicKernel::W;
	m_gradKernelFct = CubicKernel::gradW;

	m_corrs.resize(m_x.size());
	m_deltaCorrs.resize(m_x.size());
	m_densities.resize(m_x.size());
	m_factors.resize(m_x.size());
}


void SPHVolumeSampling::step(Real& avg_density_err)
{
	avg_density_err = 0.0;

	//////////////////////////////////////////////////////////////////////////
	// SPH optimization of sampling points
	//////////////////////////////////////////////////////////////////////////
	
	const int numParticles = (int)m_x.size();	
	Real V0 = m_diameter * m_diameter * m_diameter;

	const Real diameter2 = m_diameter* m_diameter;
	const Real supportRadius = m_radius* static_cast<Real>(4.0);

	START_TIMING("neighborhoodSearch");
	m_neighborhoodSearch->find_neighbors();
	STOP_TIMING_AVG

	computeDensities(m_densities, m_mass);

	// Compute cohesion => Becker 2007
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r &xi = m_x[i];
			Vector3r &corr = m_corrs[i];
			corr.setZero();

 			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
 			{
 				const unsigned int neighborIndex = getNeighbor(0, i, j);
 				const Vector3r &xj = m_x[neighborIndex];
 
 				const Vector3r xixj = xi - xj;
 				const Real r2 = xixj.dot(xixj);
				if (r2 > 1e-6)
				{
					const Real C = (xi - xj).norm() - m_diameter;
					Vector3r xixj_n = xixj;
					xixj_n.normalize();
					corr -= m_cohesion * m_mass/ m_densities[neighborIndex] * C * xixj_n * m_kernelFct(xi - xj);
				}
 			}

			if (m_adhesion != 0.0)
			{
				// Collision detection
				Vector3r cp, normal;
				const double dist = distance(xi, 0.0, normal, cp);
				if (dist < supportRadius)
				{
					corr -= m_adhesion * m_mass / m_densities[i] * (xi - cp) * m_kernelFct(xi - cp);
				}
			}
		}
	}
	


	computeDFSPHFactor();

	pressureSolveIteration(avg_density_err);

	//////////////////////////////////////////////////////////////////////////
	// Update density error
	//////////////////////////////////////////////////////////////////////////
	avg_density_err = 0.0;
	Real maxCorr = 0.0;
	for (int i = 0; i < numParticles; i++)
	{
		avg_density_err += m_densities[i] - m_density0;
		const Real corrMag = m_corrs[i].squaredNorm();
		if (corrMag > maxCorr)
			maxCorr = corrMag;
	}
	Real cfl = m_cflFactor * static_cast<Real>(0.4) * (m_diameter / (sqrt(maxCorr)));
	cfl = min(cfl, static_cast<Real>(1.0));
	avg_density_err = avg_density_err / numParticles;

	// Apply position changes
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			Vector3r &xi = m_x[i];
			xi += cfl * m_corrs[i];

			// Collision detection
			Vector3r cp, normal;
			const double dist = distance(xi, 0.0, normal, cp);
			if (dist < m_radius)
			{
				xi = cp + m_radius * normal;
			}
		}
	}
}

