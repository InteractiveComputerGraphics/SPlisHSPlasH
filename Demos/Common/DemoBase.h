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
		Shader m_shader_vector;
		Shader m_shader_scalar;
		Shader m_shader_vector_map;
		Shader m_shader_scalar_map;
		Shader m_meshShader;
		GLuint m_textureMap;
		int m_renderWalls;
		int m_colorField;
		bool m_doPause;
		Real m_pauseAt;
		Real m_stopAt;
		bool m_enablePartioExport;
		unsigned int m_framesPerSecond;
		Real m_renderMaxValue;
		Real m_renderMinValue;
		Vector3r m_oldMousePos;
		std::vector<std::vector<unsigned int>> m_selectedParticles;
		std::unique_ptr<Utilities::SceneLoader> m_sceneLoader;
		Real m_nextFrameTime;
		unsigned int m_frameCounter;
		int m_colorMapType;
		float const* m_colorMapBuffer;
		unsigned int m_colorMapLength;
#ifdef DL_OUTPUT
		Real m_nextTiming;
#endif

		virtual void initParameters();

		void initShaders();
		void initFluidData();
		void createFluidBlocks(std::map<std::string, unsigned int> &fluidIDs, std::vector<std::vector<Vector3r>> &fluidParticles, std::vector<std::vector<Vector3r>> &fluidVelocities);
		void createEmitters();

		static void selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end, void *clientData);
		static void mouseMove(int x, int y, void *clientData);
		void particleInfo();

	public:
		static int PAUSE;
		static int PAUSE_AT;
		static int STOP_AT;
		static int NUM_STEPS_PER_RENDER;
		static int PARTIO_EXPORT;
		static int PARTIO_EXPORT_FPS;
		static int RENDER_MIN_VALUE;
		static int RENDER_MAX_VALUE;
		static int RENDER_WALLS;
		static int RENDER_COLOR_FIELD;
		static int RENDER_COLOR_MAP_TYPE;

		static int ENUM_WALLS_NONE;
		static int ENUM_WALLS_PARTICLES_ALL;
		static int ENUM_WALLS_PARTICLES_NO_WALLS;
		static int ENUM_WALLS_GEOMETRY_ALL;
		static int ENUM_WALLS_GEOMETRY_NO_WALLS;

		static int ENUM_RENDER_NONE;
		static int ENUM_RENDER_VELOCITY;
		static int ENUM_RENDER_ANGULAR_VELOCITY;
		static int ENUM_RENDER_DENSITY;

		static int ENUM_COLORMAP_NONE;
		static int ENUM_COLORMAP_JET;
		static int ENUM_COLORMAP_PLASMA;

		DemoBase();
		virtual ~DemoBase();

		void init(int argc, char **argv, const char *demoName);
		void buildModel();
		void cleanup();

		void renderFluid(FluidModel *model, float *fluidColor);

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
		Shader& getShaderVector() { return m_shader_vector; }
		Shader& getShaderScalar() { return m_shader_scalar; }
		Shader& getMeshShader() { return m_meshShader; }
		void meshShaderBegin(const float *col);
		void meshShaderEnd();
		void pointShaderBegin(Shader *shader, const float *col, const bool useTexture = false);
		void pointShaderEnd(Shader *shader, const bool useTexture = false);
		Utilities::SceneLoader::Scene& getScene() { return m_scene; }

		std::vector<std::vector<unsigned int>>& getSelectedParticles() { return m_selectedParticles; }
		bool getUseParticleCaching() const { return m_useParticleCaching; }
		void setUseParticleCaching(bool val) { m_useParticleCaching = val; }

		int getColorMapType() const { return m_colorMapType; }
		void setColorMapType(const int v);
	};
}
 
#endif