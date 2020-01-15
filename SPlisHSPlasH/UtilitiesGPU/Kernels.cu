#include "Kernels.cuh"
#include "../Simulation.h"

using namespace SPH;

//////////////////////////////////////////////////////////////////
// Helper host methods 
//////////////////////////////////////////////////////////////////

KernelData::KernelData()
{
	CudaHelper::CudaMalloc(&d_W, PRECOMPUTED_KERNEL_SIZE);
	CudaHelper::CudaMalloc(&d_gradW, PRECOMPUTED_KERNEL_SIZE + 1);
}

KernelData::~KernelData()
{
	CudaHelper::CudaFree(d_W);
	CudaHelper::CudaFree(d_gradW);
}

void updateKernelData(KernelData &data)
{
	data.radius = PrecomputedKernel<CubicKernel, PRECOMPUTED_KERNEL_SIZE>::getRadius();
	data.invStepSize = PrecomputedKernel<CubicKernel, PRECOMPUTED_KERNEL_SIZE>::getInvStepSize();
	data.radius2 = data.radius * data.radius;

	CudaHelper::MemcpyHostToDevice(PrecomputedKernel<CubicKernel, PRECOMPUTED_KERNEL_SIZE>::getWeightField(), data.d_W, PRECOMPUTED_KERNEL_SIZE);
	CudaHelper::MemcpyHostToDevice(PrecomputedKernel<CubicKernel, PRECOMPUTED_KERNEL_SIZE>::getGradField(), data.d_gradW, PRECOMPUTED_KERNEL_SIZE + 1);
}

//////////////////////////////////////////////////////////////////
//Kernels for all methods 
//////////////////////////////////////////////////////////////////

__device__
Real kernelWeightPrecomputed(const Vector3r &r, const KernelData* const data)
{
	Real res = 0.0;
	const Real r2 = r.squaredNorm();
	if (r2 <= data->radius2)
	{
		const Real rl = sqrt(r2);
		//const unsigned int pos = std::min<unsigned int>((unsigned int)(rl * data->invStepSize), PRECOMPUTED_KERNEL_SIZE-2u);
		unsigned int pos = 0;
		if(static_cast<unsigned int>(rl * data->invStepSize) < PRECOMPUTED_KERNEL_SIZE-2u)
			pos = static_cast<unsigned int>(rl * data->invStepSize);
		else
			pos = PRECOMPUTED_KERNEL_SIZE-2u;
		res = static_cast<Real>(0.5)*(data->d_W[pos] + data->d_W[pos+1]);
	}
	return res;
}

__device__
Vector3r gradKernelWeightPrecomputed(const Vector3r &r, const KernelData* const data)
{
	Vector3r res;
	const Real rl = r.norm(); // rl / radius = > 0 - 1, texturSpeicher
	if (rl <= data->radius)
	{
		//const Real rl = sqrt(r2);
		//const unsigned int pos = static_cast<unsigned int>(fminf(static_cast<unsigned int>(rl * data->invStepSize), PRECOMPUTED_KERNEL_SIZE-1u));
		unsigned int pos = 0;
		if(static_cast<unsigned int>(rl * data->invStepSize) < PRECOMPUTED_KERNEL_SIZE-1u)
			pos = static_cast<unsigned int>(rl * data->invStepSize);
		else
			pos = PRECOMPUTED_KERNEL_SIZE-1u;
		res = 0.5*(data->d_gradW[pos] + data->d_gradW[pos + 1]) * r; // ersetzbar
	}
	else
		res.setZero();

	return res;
}

__device__
Real kernelWeight(const Vector3r& rin, const Real m_radius)
{
	const Real r = sqrt(rin[0] * rin[0] + rin[1] * rin[1] + rin[2] * rin[2]);
	const Real pi = 3.14159265358979323846;

	const Real h3 = m_radius*m_radius*m_radius;
	Real m_k = static_cast<Real>(8.0) / (pi*h3);
	Real m_l = static_cast<Real>(48.0) / (pi*h3);

	Real res = 0.0;
	const Real q = r / m_radius;

	if (q <= 1.0)
	{
		if (q <= 0.5)
		{
			const Real q2 = q*q;
			const Real q3 = q2*q;
			res = m_k * (static_cast<Real>(6.0)*q3 - static_cast<Real>(6.0)*q2 + static_cast<Real>(1.0));
		}
		else
		{
			res = m_k * (static_cast<Real>(2.0)*pow(static_cast<Real>(1.0) - q, 3));
		}
	}
	return res;
}

__device__
Vector3r gradKernelWeight(const Vector3r &rin, const Real m_radius)
{

	const Real pi = 3.14159265358979323846;
	const Real h3 = m_radius*m_radius*m_radius;
	const Real m_l = static_cast<Real>(48.0) / (pi*h3);

	Vector3r res;
	const Real rl = sqrt(rin[0] * rin[0] + rin[1] * rin[1] + rin[2] * rin[2]);
	const Real q = rl / m_radius;
	if ((rl > 1.0e-6) && (q <= 1.0))
	{
		const Vector3r gradq = rin * (static_cast<Real>(1.0) / (rl*m_radius));
		if (q <= 0.5)
		{
			res = m_l*q*((Real) 3.0*q - static_cast<Real>(2.0))*gradq;
		}
		else
		{
			const Real factor = static_cast<Real>(1.0) - q;
			res = m_l*(-factor*factor)*gradq;
		}
	}
	else
		res.setZero();

	return res;
}


__device__
void addForce(const Vector3r &pos, const Vector3r &f, /* output */ Vector3r* const forcesPerThread, /* output */ Vector3r* const torquesPerThread, 
	const Vector3r* const rigidBodyPositions, const uint* const forcesPerThreadIndices, const uint* const torquesPerThreadIndices, const uint index, const int id)
{
	#ifdef _OPENMP
	int tid = id;
	#else
	int tid = 0;
	#endif
	forcesPerThread[forcesPerThreadIndices[index] + tid] += f;
	torquesPerThread[torquesPerThreadIndices[index] + tid] += (pos - rigidBodyPositions[index]).cross(f);
}


__global__
void computeDensitiesGPU(/*out*/ Real* const densities, const Real* const volumes, const Real* const boundaryVolumes, const uint* const boundaryVolumeIndices, 
	const uint* const fmIndices, const Real* const densities0, const Real W_zero, const KernelData* const kernelData, 
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
 	// Boundary: Akinci2012
	const uint i = blockIdx.x * blockDim.x + threadIdx.x;

	if(i >= numParticles)
		return;

	extern __shared__ Real densities_tmp[];

	Real &density = densities_tmp[threadIdx.x];

	density = volumes[fluidModelIndex] * W_zero;
	const double3 &xi = particles[fluidModelIndex][i];

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighborsGPU(
		density += volumes[pid] * kernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
	)
	

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
  forall_boundary_neighborsGPU(
		density += boundaryVolumes[boundaryVolumeIndices[pid - nFluids] + neighborIndex] *  kernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
	)

	density *= densities0[fluidModelIndex];

	densities[fmIndices[fluidModelIndex] + i] = densities_tmp[threadIdx.x];
} 


//////////////////////////////////////////////////////////////////
//Kernels for the WCPSH method 
//////////////////////////////////////////////////////////////////

__global__
void clearAccelerationsGPU(Real* masses, Vector3r* accelerations, const Vector3r grav, const uint numActiveParticles)
{
 	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if(i >= numActiveParticles)
		return;

	// Clear accelerations of dynamic particles
	if (masses[i] != 0.0)
	{
		Vector3r &a = accelerations[i];
		a = grav;
	}
}

__global__
void updatePressureGPU(Real* const densities, const uint* const fmIndices, Real* const pressures, const Real* const densities0, const Real m_stiffness, const Real m_exponent,
	const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if(i >= numParticles)
		return;
	
	Real &density = densities[fmIndices[fluidModelIndex] + i];
	density = max(density, densities0[fluidModelIndex]);
	pressures[fmIndices[fluidModelIndex] + i] = m_stiffness * (pow(density / densities0[fluidModelIndex], m_exponent) - static_cast<Real>(1.0));
}

__global__
void computePressureAccelsGPU( /* output */ Vector3r* const pressureAccels, /* output */ Vector3r* const forcesPerThread, /* output */ Vector3r* const torquesPerThread, const uint* const forcesPerThreadIndices, 
	const uint* const torquesPerThreadIndices, const Real* const densities, const Real* const densities0, const uint* const fmIndices, const Real* const pressures, const Real* const masses, 
	const Vector3r* const rigidBodyPositions, const Real* const volumes, const Real* const boundaryVolumes, const uint* const boundaryVolumeIndices, const bool* const isDynamic, const int tid, const KernelData* kernelData,
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
   const uint i = blockIdx.x*blockDim.x + threadIdx.x;

	if(i >= numParticles)
		return;

	extern __shared__ Vector3r pressureAccels_tmp[];

	const double3 &xi = particles[fluidModelIndex][i];

	const Real density_i = densities[fmIndices[fluidModelIndex] + i];

	pressureAccels_tmp[threadIdx.x] = Vector3r(0, 0, 0);
	Vector3r &ai = pressureAccels_tmp[threadIdx.x];

	const Real dpi = pressures[fmIndices[fluidModelIndex] + i] / (density_i*density_i);
	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighborsGPU(
		const Real density_j = densities[fmIndices[pid] + neighborIndex] * densities0[fluidModelIndex] / densities0[pid];
		const Real dpj = pressures[fmIndices[pid] + neighborIndex] / (density_j*density_j);
		ai -= densities0[fluidModelIndex] * volumes[pid] * (dpi + dpj) * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	const Real dpj = pressures[fmIndices[fluidModelIndex] + i] / (densities0[fluidModelIndex] * densities0[fluidModelIndex]);
	forall_boundary_neighborsGPU(
		const Vector3r a = densities0[fluidModelIndex] * boundaryVolumes[fmIndices[pid - nFluids] + neighborIndex] * (dpi + dpj) * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
		ai -= a;
		if(isDynamic[pid - nFluids])
		{
			addForce(Vector3r(xj.x, xj.y, xj.z), masses[i] * a, forcesPerThread, torquesPerThread, rigidBodyPositions, forcesPerThreadIndices, torquesPerThreadIndices, pid - nFluids, tid);
		}
	)

	pressureAccels[i] = pressureAccels_tmp[threadIdx.x];
}

__global__ 
void updatePosPressureAccelPressureAccel(Vector3r* const positions, Vector3r* const velocities, Vector3r* const accelerations,
	const Vector3r* const pressureAccels, const Real h, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if(i >= numParticles)
		return;
	
	accelerations[i] += pressureAccels[i];
	velocities[i] += accelerations[i] * h;
	positions[i] += velocities[i] * h;
	
}


//////////////////////////////////////////////////////////////////
//Kernels for the DFSPH method 
//////////////////////////////////////////////////////////////////

__global__ 
void computeDFSPHFactors(/* out */ Real* factors, const Real* const boundaryVolumes, const uint* const boundaryVolumeIndices, const KernelData* const kernelData, 
	const unsigned int* fmIndices, const Real* fmVolumes, const Real eps,
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;
	
	Real &factor = factors[fmIndices[fluidModelIndex] + i];
	factor = 0.0;

	//////////////////////////////////////////////////////////////////////////
	// Compute gradient dp_i/dx_j * (1/k)  and dp_j/dx_j * (1/k)
	//////////////////////////////////////////////////////////////////////////

	const double3 xi = particles[fluidModelIndex][i];
	Real sum_grad_p_k = 0.0;
	Vector3r grad_p_i;
	grad_p_i.setZero();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
forall_fluid_neighborsGPU(
	const Vector3r grad_p_j = -fmVolumes[fluidModelIndex] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
	sum_grad_p_k += grad_p_j.squaredNorm();
	grad_p_i -= grad_p_j;
)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	forall_boundary_neighborsGPU(
		const Vector3r grad_p_j = -boundaryVolumes[boundaryVolumeIndices[pid - nFluids] + neighborIndex] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
		grad_p_i -= grad_p_j;
	)

	sum_grad_p_k += grad_p_i.squaredNorm();

	//////////////////////////////////////////////////////////////////////////
	// Compute pressure stiffness denominator
	//////////////////////////////////////////////////////////////////////////
	if (sum_grad_p_k > eps)
		factor = -static_cast<Real>(1.0) / (sum_grad_p_k);
	else
		factor = 0.0;
}


 __global__
void computeDensityChanges(/*out*/ Real* const densitiesAdv, const Vector3r* const fmVelocities, const Vector3r* const bmVelocities, const uint* const fmIndices, 
	const Real* const fmVolumes, const Real* const boundaryVolumes, const uint* const boundaryVolumeIndices, const KernelData* const kernelData,
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	Real &densityAdv = densitiesAdv[fmIndices[fluidModelIndex] + i];	
	const double3 &xi = particles[fluidModelIndex][i];
	const Vector3r &vi = fmVelocities[fmIndices[fluidModelIndex] + i];

	densityAdv = 0.0;
	unsigned int numNeighbors = 0;

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighborsGPU(
		const Vector3r &vj = fmVelocities[fmIndices[pid] + neighborIndex];
		densityAdv += fmVolumes[pid] * (vi - vj).dot(gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData));
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	forall_boundary_neighborsGPU(
		const Vector3r &vj = bmVelocities[boundaryVolumeIndices[pid - nFluids] + neighborIndex];
		densityAdv += boundaryVolumes[boundaryVolumeIndices[pid - nFluids] + neighborIndex] * (vi - vj).dot(gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData));
	)
	
	// only correct positive divergence
	densityAdv = max(densityAdv, static_cast<Real>(0.0));

	for (unsigned int pid = 0; pid < nPointSets; pid++)
	{
		const uint neighborsetIndex = neighborPointsetIndices[fluidModelIndex] + pid;
		numNeighbors += neighborCounts[neighborsetIndex][i];
	}

	// in case of particle deficiency do not perform a divergence solve
	if (numNeighbors < 20)
		densityAdv = 0.0;
}

__global__
void computeDensityAdvs(/*out*/ Real* const densitiesAdv, const Real* const fmDensities, const Vector3r* const fmVelocities, const Vector3r* const bmVelocities, const uint* const fmIndices, 
	const Real* const fmVolumes, const Real* const boundaryVolumes, const uint* const boundaryVolumeIndices, const Real* const densities0, const Real h, const KernelData* const kernelData,
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	Real &densityAdv = densitiesAdv[fmIndices[fluidModelIndex] + i];
	const Real &density = fmDensities[fmIndices[fluidModelIndex] + i];
	const double3 &xi = particles[fluidModelIndex][i];
	const Vector3r &vi = fmVelocities[fmIndices[fluidModelIndex] + i];
	Real delta = 0.0;

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighborsGPU(
		const Vector3r &vj = fmVelocities[fmIndices[pid] + neighborIndex];
		delta += fmVolumes[pid] * (vi - vj).dot(gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData));
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	forall_boundary_neighborsGPU(
		const Vector3r &vj = bmVelocities[boundaryVolumeIndices[pid - nFluids] + neighborIndex];
		delta += boundaryVolumes[boundaryVolumeIndices[pid - nFluids] + neighborIndex] * (vi - vj).dot(gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData));
	)
	
	densityAdv = density / densities0[fluidModelIndex] + h*delta;
	densityAdv = max(densityAdv, static_cast<Real>(1.0));
}

__global__
void warmstartDivergenceSolveKappaV(/*out*/ Real* const kappaV, const uint* const fmIndices, const Real* const densities0, const Real invH, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;
	
	kappaV[fmIndices[fluidModelIndex] + i] = static_cast<Real>(0.5) * max( kappaV[fmIndices[fluidModelIndex] + i] * invH, -static_cast<Real>(0.5) * densities0[fluidModelIndex] * densities0[fluidModelIndex]);
}

__global__
void divergenceSolveWarmstart( /*out*/ Vector3r* const fmVelocities, /* output */ Vector3r* const forcesPerThread, /* output */ Vector3r* const torquesPerThread, 
	const uint* const forcesPerThreadIndices, const uint* const torquesPerThreadIndices, const Vector3r* const rigidBodyPositions, const Real* const kappaV,
	const uint* const fmIndices, const Real* const masses, const Real* const fmVolumes, const Real* const boundaryVolumes, const uint* const boundaryVolumeIndices, 
	const Real* const densities0, const bool* const isDynamic, const int tid, const Real h, const KernelData* const kernelData, const Real eps,
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles || numParticles == 0)
		return;

	const Real invH = static_cast<Real>(1.0) / h;

	Vector3r &vel = fmVelocities[fmIndices[fluidModelIndex] + i];
	const Real ki = kappaV[fmIndices[fluidModelIndex] + i];
	const double3 &xi = particles[fluidModelIndex][i];

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighborsGPU(
		const Real kj = kappaV[fmIndices[pid] + neighborIndex];

		const Real kSum = (ki + densities0[pid] / densities0[fluidModelIndex] * kj);
		if (fabsf(kSum) > eps)
		{
			const Vector3r grad_p_j = -fmVolumes[pid] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
			vel -= h * kSum * grad_p_j;					// ki, kj already contain inverse density
		}
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (fabsf(ki) > eps)
	{
		forall_boundary_neighborsGPU(
			const Vector3r grad_p_j = -boundaryVolumes[boundaryVolumeIndices[pid - nFluids] + neighborIndex] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
			const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
			vel += velChange;
			addForce(Vector3r(xj.x, xj.y, xj.z), -masses[fmIndices[fluidModelIndex] + i] * velChange * invH, forcesPerThread, torquesPerThread, rigidBodyPositions, forcesPerThreadIndices, torquesPerThreadIndices, pid - nFluids, tid);
		)
	}
}


__global__
void multiplyRealWithConstant(/*out*/ Real* const input, const uint* const fmIndices, const Real f, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	input[fmIndices[fluidModelIndex] + i] *= f;
}

__global__
void setRealToZero(/*out*/ Real* const input, const uint* const fmIndices, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	input[fmIndices[fluidModelIndex] + i] = 0.0;
}

__global__
void divergenceSolveUpdateKappaV(/*out*/ Real* const kappaV, const Real* const densitiesAdv, const Real* const factors, const uint* const fmIndices, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	const Real b_i = densitiesAdv[fmIndices[fluidModelIndex] + i];
	const Real ki = b_i * factors[fmIndices[fluidModelIndex] + i];
	kappaV[fmIndices[fluidModelIndex] + i] += ki;
}

__global__
void pressureSolveUpdateKappa(/*out*/ Real* const kappa, const Real* const densitiesAdv, const Real* const factors, const uint* const fmIndices, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	const Real b_i = densitiesAdv[fmIndices[fluidModelIndex] + i]  - static_cast<Real>(1.0);
	const Real ki = b_i * factors[fmIndices[fluidModelIndex] + i];
	kappa[fmIndices[fluidModelIndex] + i] += ki;
}

__global__ 
void updateFluidVelocities( /*out*/ Vector3r* const fmVelocities, /* output */ Vector3r* const forcesPerThread, /* output */ Vector3r* const torquesPerThread, 
	const uint* const forcesPerThreadIndices, const uint* const torquesPerThreadIndices, const Vector3r* const rigidBodyPositions, const Real* const densitiesAdv, const Real* const factors, 
	const uint* const fmIndices, const Real* const masses, const Real* const fmVolumes, const Real* const boundaryVolumes, const uint* const boundaryVolumeIndices, 
	const Real* const densities0, const bool* const isDynamic, const int tid, const Real h, const Real invH, const KernelData* const kernelData, const Real eps,
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
 	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	const Real b_i = densitiesAdv[fmIndices[fluidModelIndex] + i];
	const Real ki = b_i * factors[fmIndices[fluidModelIndex] + i];

	Vector3r &v_i = fmVelocities[fmIndices[fluidModelIndex] + i];
	const double3 &xi = particles[fluidModelIndex][i];

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighborsGPU(
		const Real b_j = densitiesAdv[fmIndices[pid] + neighborIndex];
		const Real kj = b_j * factors[fmIndices[pid] + neighborIndex];

		const Real kSum = ki + densities0[pid] / densities0[fluidModelIndex] * kj;
		if(fabsf(kSum) > eps)
		{
			const Vector3r grad_p_j = -fmVolumes[pid] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
			v_i -= h * kSum * grad_p_j; // ki, kj already contain inverse density
		}
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if(fabsf(ki) > eps)
	{
		forall_boundary_neighborsGPU(
			const Vector3r grad_p_j = -boundaryVolumes[boundaryVolumeIndices[pid - nFluids] + neighborIndex] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
			const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;	// kj already contains inverse density
			v_i += velChange;
			addForce(Vector3r(xj.x, xj.y, xj.z), -masses[fmIndices[fluidModelIndex] + i] * velChange * invH, forcesPerThread, torquesPerThread, rigidBodyPositions, forcesPerThreadIndices, torquesPerThreadIndices, pid - nFluids, tid);
		)
	}
} 

__global__ 
void pressureSolveUpdateFluidVelocities( /*out*/ Vector3r* const fmVelocities, /* output */ Vector3r* const forcesPerThread, /* output */ Vector3r* const torquesPerThread, 
	const uint* const forcesPerThreadIndices, const uint* const torquesPerThreadIndices, const Vector3r* const rigidBodyPositions, const Real* const densitiesAdv, const Real* const factors, 
	const uint* const fmIndices, const Real* const masses, const Real* const fmVolumes, const Real* const boundaryVolumes, const uint* const boundaryVolumeIndices, 
	const Real* const densities0, const bool* const isDynamic, const int tid, const Real h, const Real invH, const KernelData* const kernelData, const Real eps,
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
 	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	const Real b_i = densitiesAdv[fmIndices[fluidModelIndex] + i] - static_cast<Real>(1.0);
	const Real ki = b_i * factors[fmIndices[fluidModelIndex] + i];

	Vector3r &v_i = fmVelocities[fmIndices[fluidModelIndex] + i];
	const double3 &xi = particles[fluidModelIndex][i];

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighborsGPU(
		const Real b_j = densitiesAdv[fmIndices[pid] + neighborIndex] - static_cast<Real>(1.0);
		const Real kj = b_j * factors[fmIndices[pid] + neighborIndex];

		const Real kSum = ki + densities0[pid] / densities0[fluidModelIndex] * kj;
		if(fabsf(kSum) > eps)
		{
			const Vector3r grad_p_j = -fmVolumes[pid] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
			v_i -= h * kSum * grad_p_j; // ki, kj already contain inverse density
		}
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if(fabsf(ki) > eps)
	{
		forall_boundary_neighborsGPU(
			const Vector3r grad_p_j = -boundaryVolumes[boundaryVolumeIndices[pid - nFluids] + neighborIndex] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
			const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;	// kj already contains inverse density
			v_i += velChange;
			addForce(Vector3r(xj.x, xj.y, xj.z), -masses[fmIndices[fluidModelIndex] + i] * velChange * invH, forcesPerThread, torquesPerThread, rigidBodyPositions, forcesPerThreadIndices, torquesPerThreadIndices, pid - nFluids, tid);
		)
	}
} 

__global__
void updateDensityErrorDivergence(/*out*/ Real* const density_errors, const Real* const densitiesAdv, const Real* const densities0, const uint* const fmIndices,
	const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	//density_errors[fluidModelIndex] += densities0[fluidModelIndex] * densitiesAdv[fmIndices[fluidModelIndex] + i];
	density_errors[0] += densities0[fluidModelIndex] * densitiesAdv[fmIndices[fluidModelIndex] + i];
}

__global__
void warmstartPressureSolveKappa(/*out*/ Real* kappa, const uint* const fmIndices, const Real* const densities0, const Real invH2, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;
	
	kappa[fmIndices[fluidModelIndex] + i] = max( kappa[fmIndices[fluidModelIndex] + i] * invH2, -static_cast<Real>(0.5) * densities0[fluidModelIndex] * densities0[fluidModelIndex]);
}

__global__
void pressureSolveWarmstart(/*out*/ Vector3r* const fmVelocities , /* output */ Vector3r* const forcesPerThread, /* output */ Vector3r* const torquesPerThread, 
	const uint* const forcesPerThreadIndices, const uint* const torquesPerThreadIndices, const Vector3r* const rigidBodyPositions,const Real* const kappa, 
	const Real* const densitiesAdv, const Real* const masses, const Real* const fmVolumes, const uint* const fmIndices, const Real* const boundaryVolumes, 
	const uint* const boundaryVolumeIndices, const Real* const densities0, const bool* const isDynamic, const int tid, const Real h, const Real eps, const KernelData* const kernelData,
	/*start of forall-parameters*/ double3** particles, uint** neighbors, uint** neighborCounts, uint** neighborOffsets, 
  uint* neighborPointsetIndices, const uint nFluids, const uint nPointSets, const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	if(densitiesAdv[fmIndices[fluidModelIndex] + i] > densities0[fluidModelIndex])
	{
		const Real invH = static_cast<Real>(1.0) / h;

		Vector3r &vel = fmVelocities[fmIndices[fluidModelIndex] + i];
		const Real &ki = kappa[fmIndices[fluidModelIndex] + i];
		const double3 &xi = particles[fluidModelIndex][i];

		//////////////////////////////////////////////////////////////////////////
		// Fluid
		//////////////////////////////////////////////////////////////////////////
		forall_fluid_neighborsGPU(
			const Real kj = kappa[fmIndices[pid] + neighborIndex];

			const Real kSum = (ki + densities0[pid] / densities0[fluidModelIndex] * kj);
			if (fabsf(kSum) > eps)
			{
				const Vector3r grad_p_j = -fmVolumes[pid] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
				vel -= h * kSum * grad_p_j;					// ki, kj already contain inverse density
			}
		)

		//////////////////////////////////////////////////////////////////////////
		// Boundary
		//////////////////////////////////////////////////////////////////////////
		if (fabsf(ki) > eps)
		{
			forall_boundary_neighborsGPU(
				const Vector3r grad_p_j = -boundaryVolumes[boundaryVolumeIndices[pid - nFluids] + neighborIndex] * gradKernelWeightPrecomputed(Vector3r(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z), kernelData);
				const Vector3r velChange = -h * (Real) 1.0 * ki * grad_p_j;				// kj already contains inverse density
				vel += velChange;
				addForce(Vector3r(xj.x, xj.y, xj.z), -masses[fmIndices[fluidModelIndex] + i] * velChange * invH, forcesPerThread, torquesPerThread, rigidBodyPositions, forcesPerThreadIndices, torquesPerThreadIndices, pid - nFluids, tid);
			)
		}
	}
}

__global__
void updateDensityErrorPressureSolve(/*out*/ Real* const density_error, const Real* const densitiesAdv, const Real* const densities0, const uint* const fmIndices,
	const uint fluidModelIndex, const uint numParticles)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= numParticles)
		return;

	//density_errors[fluidModelIndex] += densities0[fluidModelIndex] * densitiesAdv[fmIndices[fluidModelIndex] + i];
	density_error[0] += densities0[fluidModelIndex] * densitiesAdv[fmIndices[fluidModelIndex] + i] - densities0[fluidModelIndex];
}