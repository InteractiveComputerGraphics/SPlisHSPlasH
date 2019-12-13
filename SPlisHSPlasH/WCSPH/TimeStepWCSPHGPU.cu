#include "TimeStepWCSPHGPU.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataWCSPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"

#include "../../extern/cuNSearch/src/Ext_NeighborhoodSearch/src/PointSetImplementation.h"

using namespace SPH;
using namespace std;
using namespace GenParam;
using namespace cuNSearch;

int TimeStepWCSPHGPU::STIFFNESS = -1;
int TimeStepWCSPHGPU::EXPONENT = -1;

TimeStepWCSPHGPU::TimeStepWCSPHGPU() :
	TimeStep()
{
	m_simulationData.init();
	m_counter = 0;
	m_stiffness = 50.0;
	m_exponent = 7.0;

	CudaHelper::CudaMalloc(&d_kernelData, 1);

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); } });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}  
}

TimeStepWCSPHGPU::~TimeStepWCSPHGPU(void)
{
	CudaHelper::CudaFree(d_kernelData);
	CudaHelper::CudaFree(d_neighbors);
	CudaHelper::CudaFree(d_neighborCounts);
	CudaHelper::CudaFree(d_neighborOffsets);
	CudaHelper::CudaFree(d_neighborPointsetIndices); 

  Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("pressure");
		model->removeFieldByName("pressure acceleration");
	}  
}

void TimeStepWCSPHGPU::initParameters()
{
 	TimeStep::initParameters();

	STIFFNESS = createNumericParameter("stiffness", "Stiffness", &m_stiffness);
	setGroup(STIFFNESS, "WCSPH");
	setDescription(STIFFNESS, "Stiffness coefficient of EOS.");
	static_cast<RealParameter*>(getParameter(STIFFNESS))->setMinValue(1e-6);

	EXPONENT = createNumericParameter("exponent", "Exponent (gamma)", &m_exponent);
	setGroup(EXPONENT, "WCSPH");
	setDescription(EXPONENT, "Exponent of EOS.");
	static_cast<RealParameter*>(getParameter(EXPONENT))->setMinValue(1e-6); 
}

void TimeStepWCSPHGPU::initCUDA() // TODO: shift this into constructor or at best spot
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	std::vector<cuNSearch::PointSet> &pointSets = sim->getCurrent()->getPointSets();
	d_particles.resize(pointSets.size());
	for(int i = 0 ; i < pointSets.size() ; ++i)
	{
		d_particles[i] = CudaHelper::GetPointer(pointSets[i].getPointSetImplementation()->getParticles());
	}

	d_volumes.resize(nModels);
	d_densities0.resize(nModels);
	for(unsigned int pid = 0; pid < nModels; pid++)
	{
		FluidModel *fm = sim->getFluidModel(pid);
		d_volumes[pid] = fm->getVolume(0);
		d_densities0[pid] = fm->getDensity0();
	}
}

void TimeStepWCSPHGPU::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	const unsigned int nPointSets = sim->numberOfPointSets();
	TimeManager *tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

	if(!isInitialized)
	{
		initCUDA();
	}

	prepareData();

	// re-compute the precomputed kernel if necessary
	if( sim->getSupportRadius() != PrecomputedKernel<CubicKernel, PRECOMPUTED_KERNEL_SIZE>::getRadius() || !isInitialized)
	{
		PrecomputedKernel<CubicKernel, PRECOMPUTED_KERNEL_SIZE>::setRadius(sim->getSupportRadius());
		updateKernelData(kernelData);
		CudaHelper::MemcpyHostToDevice(&kernelData, d_kernelData, 1);
	
		isInitialized = true;
	}

	// for computeDensities and computePressureAccels
	thrust::device_vector<Real> d_boundaryVolumes;
	thrust::device_vector<unsigned int> d_boundaryVolumeIndices(sim->numberOfPointSets() - nModels);
	unsigned int sumBoundaryVolumes = 0;
	
	for(unsigned int pid = nModels; pid < sim->numberOfPointSets(); pid++)
	{
		BoundaryModel_Akinci2012 *bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
		d_boundaryVolumeIndices[pid - nModels] = sumBoundaryVolumes;
		sumBoundaryVolumes += bm_neighbor->getVolumes().size();
	
		d_boundaryVolumes.insert(d_boundaryVolumes.end(), bm_neighbor->getVolumes().begin(), bm_neighbor->getVolumes().end());
	}

	thrust::device_vector<unsigned int> d_fmIndices(nModels);
	std::vector<unsigned int> fmIndices(nModels); // helper to copy data back
	unsigned int sumParticles = 0;
	
	// for indexing
	for(unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		FluidModel *fm = sim->getFluidModelFromPointSet(fluidModelIndex);
		
		d_fmIndices[fluidModelIndex] = sumParticles;
	
		fmIndices[fluidModelIndex] = sumParticles;
		sumParticles += fm->getDensities().size();
	}

	Real *d_densities;
	CudaHelper::CudaMalloc( &d_densities, sumParticles);

	// Compute accelerations: a(t)
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		START_TIMING("Clearing accelerations");
		clearAccelerations(fluidModelIndex);
		STOP_TIMING_AVG;
		START_TIMING("Computing desities");

		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		std::vector<cuNSearch::PointSet> &pointSets = sim->getCurrent()->getPointSets();
		PointSetImplementation *impl = pointSets[fluidModelIndex].getPointSetImplementation();
		const Real density0 = model->getDensity0();
		const unsigned int numParticles = model->numActiveParticles();
		const Real W_zero = sim->W_zero();

		computeDensitiesGPU<<<impl->getNumberOfBlocks(), impl->getThreadsPerBlock(), impl->getThreadsPerBlock() * sizeof(Real)>>>(d_densities, CudaHelper::GetPointer(d_volumes), CudaHelper::GetPointer(d_boundaryVolumes), 
			CudaHelper::GetPointer(d_boundaryVolumeIndices), CudaHelper::GetPointer(d_fmIndices), CudaHelper::GetPointer(d_densities0), W_zero, sim->getSupportRadius(), 
			CudaHelper::GetPointer(d_particles), d_neighbors, d_neighborCounts, d_neighborOffsets, d_neighborPointsetIndices, nModels, 
			nPointSets, fluidModelIndex, numParticles);

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		CudaHelper::MemcpyDeviceToHost(d_densities, &(model->getDensity(0)), model->getDensities().size()); // TODO: does it work like this?

		STOP_TIMING_AVG;
	}

	sim->computeNonPressureForces();

	thrust::device_vector<Vector3r> d_rigidBodyPositions(sim->numberOfPointSets() - nModels);
	thrust::device_vector<bool> d_isDynamic(sim->numberOfPointSets() - nModels);
	thrust::device_vector<Vector3r> d_forcesPerThread; 
	thrust::device_vector<Vector3r> d_torquesPerThread;
	thrust::device_vector<unsigned int> d_forcesPerThreadIndices(sim->numberOfPointSets() - nModels);
	thrust::device_vector<unsigned int> d_torquesPerThreadIndices(sim->numberOfPointSets() - nModels);

	// for correct indexing
	int sumForcesPerThread = 0;
	int sumTorquesPerThread = 0;

 	for (unsigned int pid = nModels; pid < sim->numberOfPointSets(); pid++)
	{
		BoundaryModel_Akinci2012 *bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));

		d_forcesPerThread.insert(d_forcesPerThread.end(), bm_neighbor->getForcesPerThread().begin(), bm_neighbor->getForcesPerThread().end());
		d_torquesPerThread.insert(d_torquesPerThread.end(), bm_neighbor->getTorquesPerThread().begin(), bm_neighbor->getTorquesPerThread().end());

		d_forcesPerThreadIndices[pid - nModels] = sumForcesPerThread;
		d_torquesPerThreadIndices[pid - nModels] = sumTorquesPerThread;

		sumForcesPerThread += bm_neighbor->getForcesPerThread().size();
		sumTorquesPerThread += bm_neighbor->getTorquesPerThread().size();

		d_rigidBodyPositions[pid - nModels] = bm_neighbor->getRigidBodyPosition();
		d_isDynamic[pid - nModels] = bm_neighbor->isDynamic();
	}

	Real *d_pressures;
	CudaHelper::CudaMalloc(&d_pressures, sumParticles);

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		START_TIMING("Update pressure");
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		std::vector<cuNSearch::PointSet> &pointSets = sim->getCurrent()->getPointSets();
		PointSetImplementation *impl = pointSets[fluidModelIndex].getPointSetImplementation();
		const unsigned int numParticles = model->numActiveParticles();
		const unsigned int nPointSets = sim->numberOfPointSets();

		updatePressureGPU<<<impl->getNumberOfBlocks(), impl->getThreadsPerBlock()>>>( d_densities, CudaHelper::GetPointer(d_fmIndices), d_pressures, 
			CudaHelper::GetPointer(d_densities0), m_stiffness, m_exponent, fluidModelIndex, numParticles);

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		STOP_TIMING_AVG;

		START_TIMING("Compute pressure accels");

		// Pressure Accelerations
 		Vector3r *d_pressureAccels;
		CudaHelper::CudaMalloc(&d_pressureAccels, m_simulationData.getPressureAccels(fluidModelIndex).size());

		thrust::device_vector<Real> d_masses;	
		d_masses.insert(d_masses.end(), model->getMasses().begin(), model->getMasses().end());

		computePressureAccelsGPU<<<impl->getNumberOfBlocks(), impl->getThreadsPerBlock(), impl->getThreadsPerBlock() * sizeof(Vector3r)>>>( d_pressureAccels, CudaHelper::GetPointer(d_forcesPerThread), CudaHelper::GetPointer(d_torquesPerThread), CudaHelper::GetPointer(d_forcesPerThreadIndices), 
			CudaHelper::GetPointer(d_torquesPerThreadIndices), d_densities, CudaHelper::GetPointer(d_densities0), CudaHelper::GetPointer(d_fmIndices), 
			d_pressures, CudaHelper::GetPointer(d_masses), CudaHelper::GetPointer(d_rigidBodyPositions), CudaHelper::GetPointer(d_volumes), CudaHelper::GetPointer(d_boundaryVolumes), 
			CudaHelper::GetPointer(d_boundaryVolumeIndices), CudaHelper::GetPointer(d_isDynamic), omp_get_thread_num(), d_kernelData, CudaHelper::GetPointer(d_particles), d_neighbors, 
			d_neighborCounts, d_neighborOffsets, d_neighborPointsetIndices, nModels, nPointSets, fluidModelIndex, numParticles);

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		CudaHelper::MemcpyDeviceToHost( d_pressures + fmIndices[fluidModelIndex], &(m_simulationData.getPressures(fluidModelIndex)[0]), m_simulationData.getPressures(fluidModelIndex).size());
		CudaHelper::MemcpyDeviceToHost( d_pressureAccels, &(m_simulationData.getPressureAccels(fluidModelIndex)[0]), m_simulationData.getPressureAccels(fluidModelIndex).size());
		CudaHelper::CudaFree(d_pressureAccels);
		STOP_TIMING_AVG; 
	}

	CudaHelper::CudaFree(d_pressures);
	CudaHelper::CudaFree(d_densities);

	sumForcesPerThread = 0;
	sumTorquesPerThread = 0;

 	for (unsigned int pid = nModels; pid < sim->numberOfPointSets(); pid++)
	{
		BoundaryModel_Akinci2012 *bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));

		CudaHelper::MemcpyDeviceToHost( CudaHelper::GetPointer(d_forcesPerThread) + sumForcesPerThread, &(bm_neighbor->getForcesPerThread()[0]), bm_neighbor->getForcesPerThread().size());
		CudaHelper::MemcpyDeviceToHost( CudaHelper::GetPointer(d_torquesPerThread) + sumTorquesPerThread, &(bm_neighbor->getTorquesPerThread()[0]), bm_neighbor->getTorquesPerThread().size());

		sumForcesPerThread += bm_neighbor->getForcesPerThread().size();
		sumTorquesPerThread += bm_neighbor->getTorquesPerThread().size();
	}

	sim->updateTimeStepSize();

	START_TIMING("Update some quantities");
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static) 
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				if (model->getParticleState(i) == ParticleState::Active)
				{
					Vector3r &pos = model->getPosition(i);
					Vector3r &vel = model->getVelocity(i);
					Vector3r &accel = model->getAcceleration(i);
					accel += m_simulationData.getPressureAccel(fluidModelIndex, i);
					vel += accel * h;
					pos += vel * h;
				}
			}
		}
	}
	STOP_TIMING_AVG;

	sim->emitParticles();
	sim->animateParticles();

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}

void TimeStepWCSPHGPU::prepareData()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nPointSets = sim->numberOfPointSets();

	if(isInitialized)
	{
		CudaHelper::CudaFree(d_neighbors);
		CudaHelper::CudaFree(d_neighborCounts);
		CudaHelper::CudaFree(d_neighborOffsets);
		CudaHelper::CudaFree(d_neighborPointsetIndices); 
	}

	std::vector<cuNSearch::PointSet> &pointSets = sim->getCurrent()->getPointSets();

	CudaHelper::CudaMalloc(&d_neighborPointsetIndices, nPointSets);
	unsigned int neighborPointsetIndices_tmp[nPointSets];

	unsigned int neighborsetCount = 0;
	for(int i = 0 ; i < nPointSets ; ++i)
	{
		neighborPointsetIndices_tmp[i] = neighborsetCount;
		neighborsetCount += pointSets[i].n_neighborsets();	
	}

	CudaHelper::MemcpyHostToDevice(neighborPointsetIndices_tmp, d_neighborPointsetIndices, nPointSets);

	// flattened out the structures for efficiency
	CudaHelper::CudaMalloc(&d_neighbors, neighborsetCount);
	CudaHelper::CudaMalloc(&d_neighborCounts, neighborsetCount);
	CudaHelper::CudaMalloc(&d_neighborOffsets, neighborsetCount);

	for(int i = 0 ; i < nPointSets ; ++i)
	{
		const unsigned int nNeighborsets = pointSets[i].n_neighborsets();

		uint* neighbors_tmp[nNeighborsets];
		uint* neighborCounts_tmp[nNeighborsets];
		uint* neighborOffsets_tmp[nNeighborsets];

		for(int j = 0; j < nNeighborsets; j++)
		{
			neighbors_tmp[j] = pointSets[i].neighbor_indices(j);
			neighborCounts_tmp[j] = pointSets[i].neighbor_counts(j);
			neighborOffsets_tmp[j] = pointSets[i].neighbor_offsets(j);
		}

		CudaHelper::MemcpyHostToDevice(neighbors_tmp, d_neighbors + neighborPointsetIndices_tmp[i], nNeighborsets);
		CudaHelper::MemcpyHostToDevice(neighborCounts_tmp, d_neighborCounts + neighborPointsetIndices_tmp[i], nNeighborsets);
		CudaHelper::MemcpyHostToDevice(neighborOffsets_tmp, d_neighborOffsets + neighborPointsetIndices_tmp[i], nNeighborsets);
	}
}


void TimeStepWCSPHGPU::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}


void TimeStepWCSPHGPU::performNeighborhoodSearch()
{
	if (Simulation::getCurrent()->zSortEnabled())
	{
		if (m_counter % 500 == 0)
		{
			Simulation::getCurrent()->performNeighborhoodSearchSort();
			m_simulationData.performNeighborhoodSearchSort();
		}
		m_counter++;
	}

	Simulation::getCurrent()->performNeighborhoodSearch(); 
}

void TimeStepWCSPHGPU::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
 	m_simulationData.emittedParticles(model, startIndex);
 }

void TimeStepWCSPHGPU::resize()
{
 	m_simulationData.init();
}
