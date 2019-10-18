#include "TimeStep.h"
#include "TimeManager.h"
#include "SPHKernels.h"
#include "Utilities/Timing.h"
#include "EmitterSystem.h"
#include "Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;
using namespace std;
using namespace GenParam;

//#define USE_FD_NORMAL

int TimeStep::SOLVER_ITERATIONS = -1;
int TimeStep::MIN_ITERATIONS = -1;
int TimeStep::MAX_ITERATIONS = -1;
int TimeStep::MAX_ERROR = -1;


TimeStep::TimeStep()
{
	m_iterations = 0;
	m_minIterations = 2;
	m_maxIterations = 100;
	m_maxError = 0.01;
}

TimeStep::~TimeStep(void)
{
}

void TimeStep::init()
{
	initParameters();
}

void TimeStep::initParameters()
{
	ParameterObject::initParameters();

	SOLVER_ITERATIONS = createNumericParameter("iterations", "Iterations", &m_iterations);
	setGroup(SOLVER_ITERATIONS, "Simulation");
	setDescription(SOLVER_ITERATIONS, "Iterations required by the pressure solver.");
	getParameter(SOLVER_ITERATIONS)->setReadOnly(true);

	MIN_ITERATIONS = createNumericParameter("minIterations", "Min. iterations", &m_minIterations);
	setGroup(MIN_ITERATIONS, "Simulation");
	setDescription(MIN_ITERATIONS, "Minimal number of iterations of the pressure solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MIN_ITERATIONS))->setMinValue(0);

	MAX_ITERATIONS = createNumericParameter("maxIterations", "Max. iterations", &m_maxIterations);
	setGroup(MAX_ITERATIONS, "Simulation");
	setDescription(MAX_ITERATIONS, "Maximal number of iterations of the pressure solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS))->setMinValue(1);

	MAX_ERROR = createNumericParameter("maxError", "Max. density error(%)", &m_maxError);
	setGroup(MAX_ERROR, "Simulation");
	setDescription(MAX_ERROR, "Maximal density error (%).");
	static_cast<RealParameter*>(getParameter(MAX_ERROR))->setMinValue(1e-6);
}

void TimeStep::clearAccelerations(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int count = model->numActiveParticles();
	const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));
	for (unsigned int i=0; i < count; i++)
	{
		// Clear accelerations of dynamic particles
		if (model->getMass(i) != 0.0)
		{
			Vector3r &a = model->getAcceleration(i);
			a = grav;
		}
	}
}

void TimeStep::computeDensities(const unsigned int fluidModelIndex)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			Real &density = model->getDensity(i);

			// Compute current density for particle i
			density = model->getVolume(i) * sim->W_zero();
			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				density += fm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(				
					// Boundary: Akinci2012
					density += bm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					density += rho;
				);
			}
			else   // Bender2019
			{
				forall_volume_maps(
					density += Vj * sim->W(xi - xj);
				);
			}

			density *= density0;
		}
	}
}

void TimeStep::reset()
{
	m_iterations = 0;
}

void TimeStep::computeVolumeAndBoundaryX()
{
	START_TIMING("computeVolumeAndBoundaryX");
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const Vector3r &xi = model->getPosition(i);
				computeVolumeAndBoundaryX(fluidModelIndex, i, xi);
			}
		}
	}
	STOP_TIMING_AVG;
}

void TimeStep::approximateNormal(Discregrid::DiscreteGrid* map, const Eigen::Vector3d &x, Vector3r &n, const unsigned int dim)
{
	// approximate gradient
	double eps = 0.1*Simulation::getCurrent()->getSupportRadius();
	n.setZero();
	Eigen::Vector3d xTmp = x;
	for (unsigned int j = 0; j < dim; j++)
	{
		xTmp[j] += eps;

		double e_p, e_m;
		e_p = map->interpolate(0, xTmp);
		xTmp[j] = x[j] - eps;
		e_m = map->interpolate(0, xTmp);
		xTmp[j] = x[j];

		double res = (e_p - e_m) * (1.0 / (2.0*eps));

		n[j] = static_cast<Real>(res);
	}
}


void TimeStep::computeVolumeAndBoundaryX(const unsigned int fluidModelIndex, const unsigned int i, const Vector3r &xi)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const bool sim2D = sim->is2DSimulation();
	const Real supportRadius = sim->getSupportRadius();
	const Real particleRadius = sim->getParticleRadius();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	for (unsigned int pid = 0; pid < nBoundaries; pid++)
	{
		BoundaryModel_Bender2019 *bm = static_cast<BoundaryModel_Bender2019*>(sim->getBoundaryModel(pid));

		Vector3r &boundaryXj = bm->getBoundaryXj(fluidModelIndex, i);
		boundaryXj.setZero();
		Real &boundaryVolume = bm->getBoundaryVolume(fluidModelIndex, i);
		boundaryVolume = 0.0;

		const Vector3r &t = bm->getRigidBodyObject()->getPosition();
		const Matrix3r &R = bm->getRigidBodyObject()->getRotation();

		Eigen::Vector3d normal;
		const Eigen::Vector3d localXi = (R.transpose() * (xi - t)).cast<double>();


		std::array<unsigned int, 32> cell; 
		Eigen::Vector3d c0;
		Eigen::Matrix<double, 32, 1> N;
#ifdef USE_FD_NORMAL
		bool chk = bm->getMap()->determineShapeFunctions(0, localXi, cell, c0, N);
#else
		Eigen::Matrix<double, 32, 3> dN;
		bool chk = bm->getMap()->determineShapeFunctions(0, localXi, cell, c0, N, &dN);
#endif
		Real dist = numeric_limits<Real>::max();
		if (chk)
#ifdef USE_FD_NORMAL
			dist = static_cast<Real>(bm->getMap()->interpolate(0, localXi, cell, c0, N));
#else
			dist = static_cast<Real>(bm->getMap()->interpolate(0, localXi, cell, c0, N, &normal, &dN));
#endif

		if ((dist > 0.1*particleRadius) && (dist < supportRadius))
		{
			const Real volume = static_cast<Real>(bm->getMap()->interpolate(1, localXi, cell, c0, N));
			if ((volume > 1e-6) && (volume != numeric_limits<Real>::max()))
			{
				boundaryVolume = volume;

#ifdef USE_FD_NORMAL
				if (sim2D)
					approximateNormal(bm->getMap(), localXi, normal, 2);
				else
					approximateNormal(bm->getMap(), localXi, normal, 3);
#endif
				normal = R.cast<double>() * normal;
				const double nl = normal.norm();
				if (nl > 1.0e-6)
				{
					normal /= nl;
 					boundaryXj = (xi - dist * normal.cast<Real>());
				}
				else
				{
					boundaryVolume = 0.0;
				}
			}
			else
			{
				boundaryVolume = 0.0;
			}
		}
		else if (dist <= 0.1*particleRadius)
		{
			normal = R.cast<double>() * normal;
			const double nl = normal.norm();
			if (nl > 1.0e-6)
			{
				normal /= nl;
				// project to surface
				Real d = -dist;
				d = std::min(d, static_cast<Real>(0.25/0.005) * particleRadius * dt);		// get up in small steps
				sim->getFluidModel(fluidModelIndex)->getPosition(i) = (xi + d * normal.cast<Real>());
				// adapt velocity in normal direction
				sim->getFluidModel(fluidModelIndex)->getVelocity(i) += (0.05 - sim->getFluidModel(fluidModelIndex)->getVelocity(i).dot(normal.cast<Real>())) * normal.cast<Real>();
			}
			boundaryVolume = 0.0;
		}
		else
		{
			boundaryVolume = 0.0;
		}	
	}
}


void TimeStep::computeDensityAndGradient()
{
	START_TIMING("computeDensityAndGradient");
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const Vector3r &xi = model->getPosition(i);
				computeDensityAndGradient(fluidModelIndex, i, xi);
			}
		}
	}
	STOP_TIMING_AVG;
}

void TimeStep::computeDensityAndGradient(const unsigned int fluidModelIndex, const unsigned int i, const Vector3r &xi)
{
	Simulation *sim = Simulation::getCurrent();
	const bool sim2D = sim->is2DSimulation();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real supportRadius = sim->getSupportRadius();
	const Real particleRadius = sim->getParticleRadius();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	for (unsigned int pid = 0; pid < nBoundaries; pid++)
	{
		BoundaryModel_Koschier2017 *bm = static_cast<BoundaryModel_Koschier2017*>(sim->getBoundaryModel(pid));

		Vector3r &boundaryDensityGradient = bm->getBoundaryDensityGradient(fluidModelIndex, i);
		Real &boundaryDensity = bm->getBoundaryDensity(fluidModelIndex, i);

		const Vector3r &t = bm->getRigidBodyObject()->getPosition();
		const Matrix3r &R = bm->getRigidBodyObject()->getRotation();

		Eigen::Vector3d normal;
		const Eigen::Vector3d localXi = (R.transpose() * (xi - t)).cast<double>();

		Vector3r &boundaryXj = bm->getBoundaryXj(fluidModelIndex, i);		
 		std::array<unsigned int, 32> cell; 
 		Eigen::Vector3d c0;
 		Eigen::Matrix<double, 32, 1> N;
		Eigen::Matrix<double, 32, 3> dN;
		bool chk = bm->getMap()->determineShapeFunctions(0, localXi, cell, c0, N, &dN);
		Real dist = numeric_limits<Real>::max();
		if (chk)
#ifdef USE_FD_NORMAL
			dist = static_cast<Real>(bm->getMap()->interpolate(0, localXi, cell, c0, N));
#else
			dist = static_cast<Real>(bm->getMap()->interpolate(0, localXi, cell, c0, N, &normal, &dN));
#endif

		if ((dist > 0.1*particleRadius) && (dist < sim->getSupportRadius()))
		{
			Eigen::Vector3d gradD;
			const Real d = static_cast<Real>(bm->getMap()->interpolate(1, localXi, cell, c0, N, &gradD, &dN));
			if ((d > 1e-6) && (d != numeric_limits<Real>::max()))
			{
				boundaryDensity = d;

#ifdef USE_FD_NORMAL
				if (sim2D)
					approximateNormal(bm->getMap(), localXi, normal, 2);
				else
					approximateNormal(bm->getMap(), localXi, normal, 3);
#endif
				boundaryDensityGradient = -R * gradD.cast<Real>();
				normal = R.cast<double>() * normal;	
				const double nl = normal.norm();
				if (nl > 1.0e-6)
				{
					normal /= nl;
					// Move the boundary point below the surface. Otherwise
					// the particles are too close which causes large forces.
					boundaryXj = (xi - (dist+0.5*supportRadius) * normal.cast<Real>());
				}
				else
				{
					boundaryDensityGradient.setZero();
					boundaryDensity = 0.0;
				}
			}
			else
			{
				boundaryDensityGradient.setZero();
				boundaryDensity = 0.0;
			}
		}
		else if (dist <= 0.1*particleRadius)
		{
			normal = R.cast<double>() * normal;
			const double nl = normal.norm();
			if (nl > 1.0e-6)
			{
				normal /= nl;
				// project to surface
				Real d = -dist;
				d = std::min(d, static_cast<Real>(0.25 / 0.005) * particleRadius * dt);		// get up in small steps
				sim->getFluidModel(fluidModelIndex)->getPosition(i) = (xi + d * normal.cast<Real>());
				// adapt velocity in normal direction
				sim->getFluidModel(fluidModelIndex)->getVelocity(i) += (0.05 - sim->getFluidModel(fluidModelIndex)->getVelocity(i).dot(normal.cast<Real>())) * normal.cast<Real>();
			}
			boundaryDensityGradient.setZero();
			boundaryDensity = 0.0;
		}
		else
		{
			boundaryDensityGradient.setZero();
			boundaryDensity = 0.0;
		}
	}
}

