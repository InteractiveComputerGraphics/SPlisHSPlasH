#ifndef __TimeStepDFSPHGPU_h__
#define __TimeStepDFSPHGPU_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataDFSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/UtilitiesGPU/Kernels.cuh"

#include <thrust/device_vector.h>

#define USE_WARMSTART
#define USE_WARMSTART_V

namespace SPH
{
	class SimulationDataDFSPH;

	/** \brief This class implements the Divergence-free Smoothed Particle Hydrodynamics approach introduced
	* by Bender and Koschier \cite Bender:2015, \cite Bender2017, \cite KBST19.
	*/
	class TimeStepDFSPHGPU : public TimeStep
	{
	protected:
/* 		const unsigned int KERNEL_RESOLUTION = 10000;
		typedef PrecomputedKernel<CubicKernel, KERNEL_RESOLUTION> PrecomputedCubicKernel; TODO: why is this leading to compiler errors sometimes? */

		SimulationDataDFSPH m_simulationData;
		unsigned int m_counter;
		const Real m_eps = 1.0e-5;
		bool m_enableDivergenceSolver;
		unsigned int m_iterationsV;
		Real m_maxErrorV;
		unsigned int m_maxIterationsV;

		bool isInitialized = false;

		KernelData *d_kernelData, kernelData;

		thrust::device_vector<double3*> d_particles; // particle positions
		uint **d_neighbors;
		uint **d_neighborCounts;
		uint **d_neighborOffsets;
		uint *d_neighborPointsetIndices; // indexing the above

		thrust::device_vector<Real> d_volumes;
		thrust::device_vector<Real> d_densities0, d_fmDensities;
		thrust::device_vector<Vector3r> d_fmVelocities, d_bmVelocities;
		thrust::device_vector<Real> d_boundaryVolumes;
		thrust::device_vector<unsigned int> d_boundaryVolumeIndices;
		thrust::device_vector<unsigned int> d_fmIndices;

		thrust::device_vector<Real> d_masses;
		thrust::device_vector<Vector3r> d_rigidBodyPositions;
		thrust::device_vector<bool> d_isDynamic;
		thrust::device_vector<Vector3r> d_forcesPerThread; 
		thrust::device_vector<Vector3r> d_torquesPerThread;
		thrust::device_vector<unsigned int> d_forcesPerThreadIndices;
		thrust::device_vector<unsigned int> d_torquesPerThreadIndices;
		
		unsigned int sumParticles;
		
		Real *d_densitiesAdv;
		Real *d_factors;
		
	#ifdef USE_WARMSTART_V
		Real *d_kappaV;
	#endif
	
	#ifdef USE_WARMSTART
		Real *d_kappa;
	#endif

		void initCUDA();
		void prepareData();
		void getDataBack();

		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void divergenceSolveIterationDummy(const unsigned int fluidModelIndex, Real &avg_density_err);
		void divergenceSolve();
		void divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int index, const int numParticles, const Real h, const Real density0);
		void computeDensityChange(const unsigned int fluidModelIndex, const unsigned int index, const Real h);

#ifdef USE_WARMSTART_V
		void warmstartDivergenceSolve(const unsigned int fluidModelIndex);
#endif
#ifdef USE_WARMSTART
		void warmstartPressureSolve(const unsigned int fluidModelIndex);
#endif

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();
		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		virtual void initParameters();

	public:
		static int SOLVER_ITERATIONS_V;
		static int MAX_ITERATIONS_V;
		static int MAX_ERROR_V;
		static int USE_DIVERGENCE_SOLVER;

		TimeStepDFSPHGPU();
		virtual ~TimeStepDFSPHGPU(void);

		virtual void step();
		virtual void reset();

		virtual void resize();
	};
}

#endif