#ifndef __DemoBase_h__
#define __DemoBase_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "Visualization/Shader.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "extern/AntTweakBar/include/AntTweakBar.h"
#include "ParameterObject.h"

namespace SPH
{
	class DemoBase : public GenParam::ParameterObject
	{
	public: 
		struct SimulationMethod
		{
			short simulationMethod = 0;
			TimeStep *simulation = NULL;
			FluidModel model;
		};

	protected:
		unsigned int m_numberOfStepsPerRenderUpdate;
		std::string m_exePath;
		std::string m_dataPath;
		std::string m_outputPath;
		std::string m_sceneFile;
		bool m_useParticleCaching;
		Utilities::SceneLoader::Scene m_scene;
		GLint m_context_major_version;
		GLint m_context_minor_version;
		Shader m_shader;
		Shader m_meshShader;
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
		std::unique_ptr<Utilities::SceneLoader> m_sceneLoader;
		Real m_nextFrameTime;
		unsigned int m_frameCounter;


		virtual void initParameters();

		void initShaders();
		void initFluidData(std::vector<Vector3r> &fluidParticles, std::vector<Vector3r> &fluidVelocities);
		void createFluidBlocks(std::vector<Vector3r> &fluidParticles, std::vector<Vector3r> &fluidVelocities);

		static void selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end, void *clientData);
		static void mouseMove(int x, int y, void *clientData);

	public:
		static int PAUSE;
		static int PAUSE_AT;
		static int NUM_STEPS_PER_RENDER;
		static int PARTIO_EXPORT;
		static int PARTIO_EXPORT_FPS;
		static int RENDER_MAX_VEL;
		static int RENDER_OMEGA;
		static int RENDER_WALLS;

		static int ENUM_WALLS_NONE;
		static int ENUM_WALLS_PARTICLES_ALL;
		static int ENUM_WALLS_PARTICLES_NO_WALLS;
		static int ENUM_WALLS_GEOMETRY_ALL;
		static int ENUM_WALLS_GEOMETRY_NO_WALLS;

		DemoBase();
		virtual ~DemoBase();

		void init(int argc, char **argv, const char *demoName);
		void buildModel();
		void cleanup();

		void renderFluid();

		void readParameters();
		void partioExport();
		void step();
		void reset();

		Utilities::SceneLoader *getSceneLoader() { return m_sceneLoader.get(); }

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
		Utilities::SceneLoader::Scene& getScene() { return m_scene; }

		std::vector<unsigned int>& getSelectedParticles() { return m_selectedParticles; }
		bool getUseParticleCaching() const { return m_useParticleCaching; }
		void setUseParticleCaching(bool val) { m_useParticleCaching = val; }
	};
}
 
#endif