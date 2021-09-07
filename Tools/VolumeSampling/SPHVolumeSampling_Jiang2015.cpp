#include "SPlisHSPlasH/Common.h"
#include "SPHVolumeSampling_Jiang2015.h"

#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/Version.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "Utilities/FileSystem.h"


using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;

SPHVolumeSampling_Jiang2015::SPHVolumeSampling_Jiang2015() :
	SPHSamplingBase()
{
	m_stiffness = 10000.0;
	m_dt = 0.0001;
	m_cohesion = 20.0;
}

SPHVolumeSampling_Jiang2015::~SPHVolumeSampling_Jiang2015()
{

}

void SPHVolumeSampling_Jiang2015::generateSamples()
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

	std::cout << "\r" << "Generated samples: " << m_x.size() << "\n";

	initSPHOptimization();

	std::string exePath = FileSystem::getProgramPath();
	std::string exportPath = FileSystem::normalizePath(m_outputPath + "/SPHVolumeSampling_Jiang2015");
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

void SPHVolumeSampling_Jiang2015::initSPHOptimization()
{
	m_mass = static_cast<Real>(1.0)* m_diameter* m_diameter* m_diameter* m_density0;

	// Init neighborhood search
	const Real supportRadius = static_cast<Real>(4.0) * m_radius;
	m_neighborhoodSearch = new NeighborhoodSearch(supportRadius, false);
	m_neighborhoodSearch->set_radius(supportRadius);
	m_neighborhoodSearch->add_point_set(&m_x[0][0], m_x.size(), true, true);

	Poly6Kernel::setRadius(supportRadius);
	SpikyKernel::setRadius(supportRadius);
	CohesionKernel::setRadius(supportRadius);

	m_W_zero = Poly6Kernel::W_zero();
	m_kernelFct = Poly6Kernel::W;
	m_gradKernelFct = SpikyKernel::gradW;

	m_corrs.resize(m_x.size());
	m_deltaCorrs.resize(m_x.size());
	m_densities.resize(m_x.size());
	m_factors.resize(m_x.size());
	m_v.resize(m_x.size(), Vector3r::Zero());
	m_a.resize(m_x.size(), Vector3r::Zero());
	m_p.resize(m_x.size(), 0.0);
	m_n.resize(m_x.size(), Vector3r::Zero());
	//m_surfaceParticle.resize(m_x.size(), false);
}


void SPHVolumeSampling_Jiang2015::step(Real& avg_density_err)
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
	STOP_TIMING_AVG;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			m_a[i].setZero();
		}
	}

	computeDensities(m_densities, m_mass);
	computeNormals();
	computePressure();
	computeCohesion();

	// time integration
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			m_v[i] += m_dt * m_a[i];

			// damping 
			m_v[i] *= 0.9;

			Vector3r cp, normal;
			const double dist = distance(m_x[i], 0.0, normal, cp);
			if (dist < 0.0)
			{
				m_v[i] -= m_v[i].dot(m_n[i]) * m_n[i];
			}
			else
				m_x[i] += m_dt * m_v[i];
		}
	}
}

void SPHVolumeSampling_Jiang2015::computePressure()
{
	const int numParticles = (int)m_x.size();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			// compute pressure (EOS)
			m_p[i] = std::max(m_stiffness * (m_densities[i] - m_density0), static_cast<Real>(0.0));
		}

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			// compute pressure accel
			const Vector3r& xi = m_x[i];
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r& xj = m_x[neighborIndex];

				m_a[i] -= m_mass * (m_p[i]/(m_densities[i]*m_densities[i]) + m_p[neighborIndex] / (m_densities[neighborIndex] * m_densities[neighborIndex])) * m_gradKernelFct(xi - xj);
			}

			Vector3r cp, normal;
			const double dist = distance(xi, 0.0, normal, cp);
			if (dist <= 0.0)
			{
				// remove normal acceleration
				m_a[i] -= m_a[i].dot(m_n[i]) * m_n[i];
			}
		}

	}
}

void SPHVolumeSampling_Jiang2015::computeCohesion()
{
	const int numParticles = (int)m_x.size();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			// compute pressure accel
			const Vector3r& xi = m_x[i];
			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r& xj = m_x[neighborIndex];

				const Real Kij = static_cast<Real>(2.0) * m_density0 / (m_densities[i] + m_densities[neighborIndex]);
				Vector3r xixj = xi - xj;
				xixj.normalize();
				m_a[i] -= m_cohesion * m_mass * Kij * CohesionKernel::W(xi - xj) * xixj;
			}
		}

	}
}


void SPHVolumeSampling_Jiang2015::computeNormals()
{
	const int numParticles = (int)m_x.size();

	// Compute normals
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r& xi = m_x[i];
			Vector3r &ni = m_n[i];
			ni.setZero();

			for (unsigned int j = 0; j < numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(0, i, j);
				const Vector3r& xj = m_x[neighborIndex];

				ni -= m_mass / m_densities[neighborIndex] * m_gradKernelFct(xi - xj);
			}
			ni.normalize();
		}
	}
}
