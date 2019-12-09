#ifndef __TimeStepWCSPHGPU_h__
#define __TimeStepWCSPHGPU_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataWCSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/UtilitiesGPU/Kernels.cuh"
//#include "../Common.h"

#include <thrust/device_vector.h>

namespace SPH
{
	class SimulationDataWCSPH;

/* 	struct KernelDeviceData{
		Real *m_W;
		Real *m_gradW;
		Real m_radius;
		Real m_radius2;
		Real m_invStepSize;
		Real m_W_zero;
	}; */

	/** \brief This class implements the Weakly Compressible SPH for Free Surface Flows approach introduced
	* by Becker and Teschner \cite Becker:2007.
	*/
	class TimeStepWCSPHGPU : public TimeStep
	{
	protected:

/* 		KernelDeviceData kernelData, *d_kernelData;
		unsigned int kernelResolution; */

		bool isInitialized = false;

		KernelData *d_kernelData, kernelData;

		thrust::device_vector<double3*> d_particles; // particle positions
	  thrust::device_vector<uint> d_neighborIndices; // indices of the neighbor particles
	  thrust::device_vector<uint> d_neighborCounts; // how many neighbors
	  thrust::device_vector<uint> d_neighborWriteOffsets; // offsets, where neighbors within neid_neighborIndices begin

		thrust::device_vector<uint> d_psPidStartIndices; // where a pid starts within neihborcounts
		thrust::device_vector<uint> d_neighborPidStartIndices; // where a pid starts within neighborIndices
		thrust::device_vector<uint> d_pointsetIndices; // where a point set starts within psPidStartIndices

		thrust::device_vector<Real> d_volumes;
		thrust::device_vector<Real> d_densities0;

		Real m_stiffness;
		Real m_exponent;

		SimulationDataWCSPH m_simulationData;
		unsigned int m_counter;

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);
		virtual void initParameters();
		void initCUDA();

		void computePressureAccels(const unsigned int fluidModelIndex);

	public:
		static int STIFFNESS;
		static int EXPONENT;

		TimeStepWCSPHGPU();
		virtual ~TimeStepWCSPHGPU(void);

		virtual void step();
		virtual void reset();
		virtual void resize();
	};
}


#endif