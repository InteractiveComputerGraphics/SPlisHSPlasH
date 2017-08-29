#ifndef __DemoBase_h__
#define __DemoBase_h__

#include "SPlisHSPlasH/Common.h"
#include "Utilities/SceneLoader.h"
#include "Visualization/Shader.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "extern/AntTweakBar/include/AntTweakBar.h"

namespace SPH
{
	class DemoBase
	{
	public: 
		struct SimulationMethod
		{
			short simulationMethod = 0;
			TimeStep *simulation = NULL;
			FluidModel model;
		};

		struct Parameter
		{
			unsigned int id;
			std::string name;
			std::string tweakBarDefinition;
			TwType type;
			DemoBase *base;

			Parameter(const unsigned int pId, const std::string &pName, const TwType pType, const std::string &pTweakBarDefinition, DemoBase *pBase) :
				id(pId), name(pName), type(pType), tweakBarDefinition(pTweakBarDefinition), base(pBase) {}
		};

		enum ParameterIDs {
			Separator = 0, TimeStepSize, IterationCount, IterationCountV,
			Gravitation, SimMethod, VelocityUpdateMethod, 
			Vorticity, VorticityCoeff, ViscosityOmega, InertiaInverse, BaroclinicCoeff, 
			DragMethod, DragCoefficient,
			BuoyancyMethod, BuoyancyCoefficient, ThermalConductivity, RadiationHalfTime, AmbientAirPressure, 
			Viscosity, ViscosityMethod, 
			ViscoMaxIter, ViscoMaxError,
			Stiffness, WCSPH_Exponent,
			DFSPH_EnableDivergenceSolver,
			CFL_Method, CFL_Factor, CFL_MaxTimeStepSize, 
			Kernel_Method, GradKernel_Method, 
			SurfaceTension, SurfaceTensionMethod,
			MaxIterations, MaxError, MaxIterationsV, MaxErrorV,
			NumParticles, ReusedParticles
		};

		enum SimulationMethods { WCSPH = 0, PCISPH, PBF, IISPH, DFSPH, PF, NUM_METHODS };

		typedef void(*SimulationMethodChangedFct)();

	protected:
		unsigned int m_numberOfStepsPerRenderUpdate;
		std::string m_exePath;
		std::string m_dataPath;
		std::string m_sceneFile;
		bool m_useParticleCaching;
		SceneLoader::Scene m_scene;
		GLint m_context_major_version;
		GLint m_context_minor_version;
		Shader m_shader;
		Shader m_meshShader;
		SimulationMethod m_simulationMethod;
		std::vector<Parameter> m_parameters;
		int m_renderWalls;
		bool m_doPause;
		Real m_pauseAt;
		bool m_enablePartioExport;
		unsigned int m_framesPerSecond;
		bool m_renderAngularVelocities;
		bool m_renderTemperatures;
		Real m_renderMaxVelocity;
		Vector3r m_oldMousePos;
		std::vector<unsigned int> m_selectedParticles;
		SimulationMethodChangedFct m_simulationMethodChangedFct;

		void initShaders();
		void initParameters();
		void initFluidData(std::vector<Vector3r> &fluidParticles, std::vector<Vector3r> &fluidVelocities);
		void createFluidBlocks(std::vector<Vector3r> &fluidParticles, std::vector<Vector3r> &fluidVelocities);

		static void TW_CALL setParameter(const void *value, void *clientData);
		static void TW_CALL getParameter(void *value, void *clientData);

		static void selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end, void *clientData);
		static void mouseMove(int x, int y, void *clientData);

	public:
		DemoBase();
		virtual ~DemoBase();

		void init(int argc, char **argv, const char *demoName);
		void buildModel();
		void cleanup();

		void renderFluid();

		unsigned int getNumberOfStepsPerRenderUpdate() const { return m_numberOfStepsPerRenderUpdate; }
		void setNumberOfStepsPerRenderUpdate(unsigned int val) { m_numberOfStepsPerRenderUpdate = val; }

		const std::string& getExePath() const { return m_exePath; }
		const std::string& getDataPath() const { return m_dataPath; }
		const std::string& getSceneFile() const { return m_sceneFile; }

		GLint getContextMajorVersion() const { return m_context_major_version; }
		GLint getContextMinorVersion() const { return m_context_minor_version; }
		Shader& getShader() { return m_shader; }
		Shader& getMeshShader() { return m_meshShader; }
		void meshShaderBegin(const float *col);
		void meshShaderEnd();
		void pointShaderBegin(const float *col);
		void pointShaderEnd();
		SceneLoader::Scene& getScene() { return m_scene; }
		SimulationMethod &getSimulationMethod() { return m_simulationMethod; }

		int getRenderWalls() const { return m_renderWalls; }
		void setRenderWalls(int val) { m_renderWalls = val; }
		bool getPause() const { return m_doPause; }
		void setPause(bool val) { m_doPause = val; }
		std::vector<unsigned int>& getSelectedParticles() { return m_selectedParticles; }
		bool getUseParticleCaching() const { return m_useParticleCaching; }
		void setUseParticleCaching(bool val) { m_useParticleCaching = val; }
		SPH::DemoBase::SimulationMethodChangedFct getSimulationMethodChangedFct() const { return m_simulationMethodChangedFct; }
		void setSimulationMethodChangedFct(SPH::DemoBase::SimulationMethodChangedFct val) { m_simulationMethodChangedFct = val; }
		Real getPauseAt() const { return m_pauseAt; }
		void setPauseAt(Real val) { m_pauseAt = val; }
		void setSimulationMethod(SimulationMethods method);
		bool getEnablePartioExport() const { return m_enablePartioExport; }
		void setEnablePartioExport(bool val) { m_enablePartioExport = val; }
		unsigned int getFramesPerSecond() const { return m_framesPerSecond; }
		void setFramesPerSecond(unsigned int val) { m_framesPerSecond = val; }
		Real getRenderMaxVelocity() const { return m_renderMaxVelocity; }
		void setRenderMaxVelocity(Real val) { m_renderMaxVelocity = val; }
		bool getRenderAngularVelocities() const { return m_renderAngularVelocities; }
		void setRenderAngularVelocities(bool val) { m_renderAngularVelocities = val; }
		bool getRenderTemperatures() const { return m_renderTemperatures; }
		void setRenderTemperatures(bool val) { m_renderTemperatures = val; }
	};
}
 
#endif