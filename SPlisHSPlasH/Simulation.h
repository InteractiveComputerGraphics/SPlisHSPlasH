#ifndef __Simulation_h__
#define __Simulation_h__

#include "Common.h"
#include "FluidModel.h"
#include "NonPressureForceBase.h"
#include "ParameterObject.h"
#include "NeighborhoodSearch.h"
#include "BoundaryModel.h"
#include "AnimationFieldSystem.h"
#include "Utilities/FileSystem.h"
#ifdef USE_DEBUG_TOOLS
#include "SPlisHSPlasH/Utilities/DebugTools.h"
#endif
#include <array>
#include <algorithm>


/** Loop over the fluid neighbors of all fluid phases. 
* Simulation *sim and unsigned int fluidModelIndex must be defined.
*/
#define forall_fluid_neighbors(code) \
	for (unsigned int pid = 0; pid < nFluids; pid++) \
	{ \
		FluidModel *fm_neighbor = sim->getFluidModelFromPointSet(pid); \
		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++) \
		{ \
			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j); \
			const Vector3r &xj = fm_neighbor->getPosition(neighborIndex); \
			code \
		} \
	} 

/** Loop over the fluid neighbors of the same fluid phase.
* Simulation *sim, unsigned int fluidModelIndex and FluidModel* model must be defined.
*/
#define forall_fluid_neighbors_in_same_phase(code) \
	for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++) \
	{ \
		const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j); \
		const Vector3r &xj = model->getPosition(neighborIndex); \
		code \
	} 

/** Loop over the boundary neighbors of all fluid phases.
* Simulation *sim and unsigned int fluidModelIndex must be defined.
*/
#define forall_boundary_neighbors(code) \
for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) \
{ \
	BoundaryModel_Akinci2012 *bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid)); \
	for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++) \
	{ \
		const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j); \
		const Vector3r &xj = bm_neighbor->getPosition(neighborIndex); \
		code \
	} \
}

/** Loop over the boundary density maps.
* Simulation *sim, unsigned int nBoundaries and unsigned int fluidModelIndex must be defined.
*/
#define forall_density_maps(code) \
for (unsigned int pid = 0; pid < nBoundaries; pid++) \
{ \
	BoundaryModel_Koschier2017 *bm_neighbor = static_cast<BoundaryModel_Koschier2017*>(sim->getBoundaryModel(pid)); \
	const Real rho = bm_neighbor->getBoundaryDensity(fluidModelIndex, i); \
	if (rho != 0.0) \
	{ \
		const Vector3r &gradRho = bm_neighbor->getBoundaryDensityGradient(fluidModelIndex, i).cast<Real>(); \
		const Vector3r &xj = bm_neighbor->getBoundaryXj(fluidModelIndex, i); \
		code \
	} \
}

/** Loop over the boundary volume maps.
* Simulation *sim, unsigned int nBoundaries and unsigned int fluidModelIndex must be defined.
*/
#define forall_volume_maps(code) \
for (unsigned int pid = 0; pid < nBoundaries; pid++) \
{ \
	BoundaryModel_Bender2019 *bm_neighbor = static_cast<BoundaryModel_Bender2019*>(sim->getBoundaryModel(pid)); \
	const Real Vj = bm_neighbor->getBoundaryVolume(fluidModelIndex, i);  \
	if (Vj > 0.0) \
	{ \
		const Vector3r &xj = bm_neighbor->getBoundaryXj(fluidModelIndex, i); \
		code \
	} \
}

#ifdef USE_AVX
/** Loop over the fluid neighbors of all fluid phases.
* Simulation *sim and unsigned int fluidModelIndex must be defined.
*/
#define forall_fluid_neighbors_avx(code) \
	for (unsigned int pid = 0; pid < nFluids; pid++) \
	{ \
		FluidModel *fm_neighbor = sim->getFluidModelFromPointSet(pid); \
		const unsigned int maxN = sim->numberOfNeighbors(fluidModelIndex, pid, i); \
		for (unsigned int j = 0; j < maxN; j += 8) \
		{ \
			const unsigned int count = std::min(maxN - j, 8u); \
			const Vector3f8 xj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getPosition(0), count); \
			code \
		} \
	} 

/** Loop over the fluid neighbors of all fluid phases.
* Simulation *sim and unsigned int fluidModelIndex must be defined.
*/
#define forall_fluid_neighbors_avx_nox(code) \
	unsigned int idx = 0; \
	for (unsigned int pid = 0; pid < nFluids; pid++) \
	{ \
		FluidModel *fm_neighbor = sim->getFluidModelFromPointSet(pid); \
		const unsigned int maxN = sim->numberOfNeighbors(fluidModelIndex, pid, i); \
		for (unsigned int j = 0; j < maxN; j += 8) \
		{ \
			const unsigned int count = std::min(maxN - j, 8u); \
			code \
			idx++; \
		} \
	} 

/** Loop over the fluid neighbors of the same fluid phase.
* Simulation *sim, unsigned int fluidModelIndex and FluidModel* model must be defined.
*/
#define forall_fluid_neighbors_in_same_phase_avx(code) \
	const unsigned int maxN = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); \
	for (unsigned int j = 0; j < maxN; j += 8) \
	{ \
		const unsigned int count = std::min(maxN - j, 8u); \
		const Vector3f8 xj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getPosition(0), count); \
		code \
	} 

/** Loop over the fluid neighbors of the same fluid phase.
* Simulation *sim, unsigned int fluidModelIndex and FluidModel* model must be defined.
*/
#define forall_fluid_neighbors_in_same_phase_avx_nox(code) \
	const unsigned int maxN = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); \
	for (unsigned int j = 0; j < maxN; j += 8) \
	{ \
		const unsigned int count = std::min(maxN - j, 8u); \
		code \
	} 

/** Loop over the boundary neighbors of all fluid phases.
* Simulation *sim and unsigned int fluidModelIndex must be defined.
*/
#define forall_boundary_neighbors_avx(code) \
	for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) \
	{ \
		BoundaryModel_Akinci2012 *bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid)); \
		const unsigned int maxN = sim->numberOfNeighbors(fluidModelIndex, pid, i); \
		for (unsigned int j = 0; j < maxN; j += 8) \
		{ \
			const unsigned int count = std::min(maxN - j, 8u); \
			const Vector3f8 xj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getPosition(0), count); \
			code \
		} \
	}

#endif



namespace SPH
{
	enum class SimulationMethods { WCSPH = 0, PCISPH, PBF, IISPH, DFSPH, PF, ICSPH, NumSimulationMethods };
	enum class BoundaryHandlingMethods { Akinci2012 = 0, Koschier2017, Bender2019, NumSimulationMethods };

	/** \brief Class to manage the current simulation time and the time step size. 
	* This class is a singleton.
	*/
	class Simulation : public GenParam::ParameterObject
	{
	public:
		static int SIM_2D;
		static int PARTICLE_RADIUS;
		static int GRAVITATION;
		static int CFL_METHOD;
		static int CFL_FACTOR;
		static int CFL_MIN_TIMESTEPSIZE;
		static int CFL_MAX_TIMESTEPSIZE;
		static int ENABLE_Z_SORT;
		static int STEPS_PER_Z_SORT;

		static int KERNEL_METHOD;
		static int GRAD_KERNEL_METHOD;
		static int ENUM_KERNEL_CUBIC;
		static int ENUM_KERNEL_WENDLANDQUINTICC2;
		static int ENUM_KERNEL_POLY6;
		static int ENUM_KERNEL_SPIKY;
		static int ENUM_KERNEL_PRECOMPUTED_CUBIC;
		static int ENUM_KERNEL_CUBIC_2D;
		static int ENUM_KERNEL_WENDLANDQUINTICC2_2D;
		static int ENUM_GRADKERNEL_CUBIC;
		static int ENUM_GRADKERNEL_WENDLANDQUINTICC2;
		static int ENUM_GRADKERNEL_POLY6;
		static int ENUM_GRADKERNEL_SPIKY;
		static int ENUM_GRADKERNEL_PRECOMPUTED_CUBIC;
		static int ENUM_GRADKERNEL_CUBIC_2D;
		static int ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D;

		static int SIMULATION_METHOD;

		static int ENUM_CFL_NONE;
		static int ENUM_CFL_STANDARD;
		static int ENUM_CFL_ITER;

		static int ENUM_SIMULATION_WCSPH;
		static int ENUM_SIMULATION_PCISPH;
		static int ENUM_SIMULATION_PBF;
		static int ENUM_SIMULATION_IISPH;
		static int ENUM_SIMULATION_DFSPH;
		static int ENUM_SIMULATION_PF;
		static int ENUM_SIMULATION_ICSPH;

		static int BOUNDARY_HANDLING_METHOD;
		static int ENUM_AKINCI2012;
		static int ENUM_KOSCHIER2017;
		static int ENUM_BENDER2019;

		typedef PrecomputedKernel<CubicKernel, 10000> PrecomputedCubicKernel;

		struct NonPressureForceMethod
		{
			std::string m_name;
			std::function<NonPressureForceBase* (FluidModel*)> m_creator;
			int m_id;
		};

		/** Fluid object information */
		struct FluidInfo
		{
			int type;		// 0: block, 1: fluid model, 2: emitter
			int numParticles;
			AlignedBox3r box;
			std::string id;
			std::string samplesFile;
			std::string visMeshFile;
			Vector3r translation;
			Matrix3r rotation;
			Vector3r scale;
			Vector3r initialVelocity;
			Vector3r initialAngularVelocity;
			unsigned char mode;
			bool invert;
			std::array<unsigned int, 3> resolutionSDF;
			unsigned int emitter_width;
			unsigned int emitter_height;
			Real emitter_velocity; // emission velocity
			Real emitter_emitStartTime;
			Real emitter_emitEndTime;
			unsigned int emitter_type;

			bool hasSameParticleSampling(const FluidInfo &other)
			{
				if (numParticles != other.numParticles)
					return false;
				if ((type == 1) && (other.type == 1) && (scale == other.scale))
				{
					if (samplesFile == other.samplesFile)
					{
						std::string ext = Utilities::FileSystem::getFileExt(samplesFile);
						std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
						if (ext == "OBJ")
						{
							if (mode == other.mode)
								return true;
						}
						else
							return true;
					}
				}
				else if ((type == 0) && (other.type == 0))
				{
					if (((box.max() - box.min()).isApprox(other.box.max() - other.box.min()), 1.0e-9) && (mode == other.mode))
						return true;
				}
				return false;
			}
		};

	protected:
		std::vector<FluidModel*> m_fluidModels;
		std::vector<BoundaryModel*> m_boundaryModels;
		std::vector<FluidInfo> m_fluidInfos;
		NeighborhoodSearch *m_neighborhoodSearch;
		AnimationFieldSystem *m_animationFieldSystem;
		int m_cflMethod;
		Real m_cflFactor;
		Real m_cflMinTimeStepSize;
		Real m_cflMaxTimeStepSize;
		int m_kernelMethod;
		int m_gradKernelMethod;
		Real m_W_zero;
		Real(*m_kernelFct)(const Vector3r &);
		Vector3r(*m_gradKernelFct)(const Vector3r &r);
		SimulationMethods m_simulationMethod;
		TimeStep *m_timeStep;
		Vector3r m_gravitation;
		Real m_particleRadius;
		Real m_supportRadius;
		bool m_sim2D;
		bool m_enableZSort;
		unsigned int m_stepsPerZSort;
		unsigned int m_counter;
		std::function<void()> m_simulationMethodChanged;		
		int m_boundaryHandlingMethod;
		std::string m_cachePath;
		bool m_useCache;
		std::vector<NonPressureForceMethod> m_dragMethods;
		std::vector<NonPressureForceMethod> m_elasticityMethods;
		std::vector<NonPressureForceMethod> m_surfaceTensionMethods;
		std::vector<NonPressureForceMethod> m_vorticityMethods;
		std::vector<NonPressureForceMethod> m_viscoMethods;
		bool m_simulationIsInitialized;
#ifdef USE_DEBUG_TOOLS
		DebugTools* m_debugTools;
#endif

		virtual void initParameters();

		void registerNonpressureForces();
		
	private:
		static Simulation *current;

	public:
		Simulation ();
		Simulation(const Simulation&) = delete;
        Simulation& operator=(const Simulation&) = delete;
		~Simulation ();

		void init(const Real particleRadius, const bool sim2D);

		/** This function is called after the simulation scene is loaded and all
		* parameters are initialized. While reading a scene file several parameters 
		* can change. The deferred init function should initialize all values which 
		* depend on these parameters. 
		*/
		void deferredInit();
		void reset();

		// Singleton
		static Simulation* getCurrent ();
		static void setCurrent (Simulation* tm);
		static bool hasCurrent();

		void addFluidModel(const std::string &id, const unsigned int nFluidParticles, Vector3r* fluidParticles, Vector3r* fluidVelocities, unsigned int* fluidObjectIds, const unsigned int nMaxEmitterParticles);
		FluidModel *getFluidModel(const unsigned int index) { return m_fluidModels[index]; }
		FluidModel *getFluidModelFromPointSet(const unsigned int pointSetIndex) { return static_cast<FluidModel*>(m_neighborhoodSearch->point_set(pointSetIndex).get_user_data()); }
		const unsigned int numberOfFluidModels() const { return static_cast<unsigned int>(m_fluidModels.size()); }

		void addBoundaryModel(BoundaryModel *bm);
		BoundaryModel *getBoundaryModel(const unsigned int index) { return m_boundaryModels[index]; }
		BoundaryModel *getBoundaryModelFromPointSet(const unsigned int pointSetIndex) { return static_cast<BoundaryModel*>(m_neighborhoodSearch->point_set(pointSetIndex).get_user_data()); }
		const unsigned int numberOfBoundaryModels() const { return static_cast<unsigned int>(m_boundaryModels.size()); }
		void updateBoundaryVolume();

		void addFluidInfo(FluidInfo& info) { m_fluidInfos.push_back(info); }
		std::vector<FluidInfo>& getFluidInfos() { return m_fluidInfos; }
		FluidInfo& getFluidInfo(const unsigned int i) { return m_fluidInfos[i]; }

		AnimationFieldSystem* getAnimationFieldSystem() { return m_animationFieldSystem; }
		
		BoundaryHandlingMethods getBoundaryHandlingMethod() const { return (BoundaryHandlingMethods) m_boundaryHandlingMethod; }
		void setBoundaryHandlingMethod(BoundaryHandlingMethods val) { m_boundaryHandlingMethod = (int) val; }

		int getKernel() const { return m_kernelMethod; }
		void setKernel(int val);
		int getGradKernel() const { return m_gradKernelMethod; }
		void setGradKernel(int val);

		int isSimulationInitialized() const { return m_simulationIsInitialized; }
		void setSimulationInitialized(int val);

		FORCE_INLINE Real W_zero() const { return m_W_zero; }
		FORCE_INLINE Real W(const Vector3r &r) const { return m_kernelFct(r); }
		FORCE_INLINE Vector3r gradW(const Vector3r &r) { return m_gradKernelFct(r); }

		int getSimulationMethod() const { return static_cast<int>(m_simulationMethod); }
		void setSimulationMethod(const int val);

		void setSimulationMethodChangedCallback(std::function<void()> const& callBackFct);

		TimeStep *getTimeStep() { return m_timeStep; }

		bool is2DSimulation() { return m_sim2D; }
		bool zSortEnabled() { return m_enableZSort; }

		unsigned int stepsPerZSort() { return m_stepsPerZSort; }

		void initKernels();

		void setParticleRadius(Real val);
		Real getParticleRadius() const { return m_particleRadius; }
		Real getSupportRadius() const { return m_supportRadius; }

		/** Update time step size depending on the chosen method.
		*/
		void updateTimeStepSize();

		/** Update time step size by CFL condition.
		*/
		void updateTimeStepSizeCFL();

		/** Perform the neighborhood search for all fluid particles.
		*/
		virtual void performNeighborhoodSearch();
		void performNeighborhoodSearchSort();

		void computeNonPressureForces();

		void animateParticles();
		void emitParticles();
		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		NeighborhoodSearch* getNeighborhoodSearch() { return m_neighborhoodSearch; }

		void setCachePath(const std::string& cachePath) { m_cachePath = cachePath; }
		const std::string& getCachePath() const { return m_cachePath; }
		void setUseCache(const bool useCache) { m_useCache = useCache; }
		const bool getUseCache() const { return m_useCache; }

		void saveState(BinaryFileWriter &binWriter);
		void loadState(BinaryFileReader &binReader);

		void addDragMethod(const std::string& name, const std::function<NonPressureForceBase* (FluidModel*)>& creator) { m_dragMethods.push_back({ name, creator, -1 }); }
		std::vector<NonPressureForceMethod>& getDragMethods() { return m_dragMethods; }

		void addElasticityMethod(const std::string& name, const std::function<NonPressureForceBase* (FluidModel*)>& creator) { m_elasticityMethods.push_back({ name, creator, -1 }); }
		std::vector<NonPressureForceMethod>& getElasticityMethods() { return m_elasticityMethods; }

		void addSurfaceTensionMethod(const std::string& name, const std::function<NonPressureForceBase* (FluidModel*)>& creator) { m_surfaceTensionMethods.push_back({ name, creator, -1 }); }
		std::vector<NonPressureForceMethod>& getSurfaceTensionMethods() { return m_surfaceTensionMethods; }

		void addViscosityMethod(const std::string& name, const std::function<NonPressureForceBase* (FluidModel*)>& creator) { m_viscoMethods.push_back({ name, creator, -1 }); }
		std::vector<NonPressureForceMethod>& getViscosityMethods() { return m_viscoMethods; }

		void addVorticityMethod(const std::string& name, const std::function<NonPressureForceBase* (FluidModel*)>& creator) { m_vorticityMethods.push_back({ name, creator, -1 }); }
		std::vector<NonPressureForceMethod>& getVorticityMethods() { return m_vorticityMethods; }

#ifdef USE_DEBUG_TOOLS
		DebugTools* getDebugTools() { return m_debugTools; }
		void createDebugTools() { m_debugTools = new DebugTools(); m_debugTools->init(); }
#endif

		FORCE_INLINE unsigned int numberOfPointSets() const
		{
			return static_cast<unsigned int>(m_neighborhoodSearch->n_point_sets());
		}

		FORCE_INLINE unsigned int numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index) const
		{
			return static_cast<unsigned int>(m_neighborhoodSearch->point_set(pointSetIndex).n_neighbors(neighborPointSetIndex, index));
		}

		FORCE_INLINE unsigned int getNeighbor(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index, const unsigned int k) const
		{
			return m_neighborhoodSearch->point_set(pointSetIndex).neighbor(neighborPointSetIndex, index, k);
		}

		FORCE_INLINE const unsigned int * getNeighborList(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index) const
		{
			#ifdef GPU_NEIGHBORHOOD_SEARCH
			return m_neighborhoodSearch->point_set(pointSetIndex).neighbor_list(neighborPointSetIndex, index);
			#else
			return m_neighborhoodSearch->point_set(pointSetIndex).neighbor_list(neighborPointSetIndex, index).data();
			#endif
		}
	};
}

#endif
