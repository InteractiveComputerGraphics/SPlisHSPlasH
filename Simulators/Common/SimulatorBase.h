#ifndef __SimulatorBase_h__
#define __SimulatorBase_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "Visualization/Shader.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "extern/AntTweakBar/include/AntTweakBar.h"
#include "ParameterObject.h"

namespace SPH
{
	class SimulatorBase : public GenParam::ParameterObject
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
		bool m_doPause;
		Real m_pauseAt;
		Real m_stopAt;
		bool m_enablePartioExport;
		unsigned int m_framesPerSecond;
		std::string m_partioAttributes;
		Vector3r m_oldMousePos;
		std::vector<std::vector<unsigned int>> m_selectedParticles;
		std::unique_ptr<Utilities::SceneLoader> m_sceneLoader;
		Real m_nextFrameTime;
		unsigned int m_frameCounter;
		std::vector<std::string> m_colorField;
		std::vector<int> m_colorMapType;
		std::vector<Real> m_renderMaxValue;
		std::vector<Real> m_renderMinValue;
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
		static int PARTIO_EXPORT_ATTRIBUTES;
		static int RENDER_WALLS;
		
		static int ENUM_WALLS_NONE;
		static int ENUM_WALLS_PARTICLES_ALL;
		static int ENUM_WALLS_PARTICLES_NO_WALLS;
		static int ENUM_WALLS_GEOMETRY_ALL;
		static int ENUM_WALLS_GEOMETRY_NO_WALLS;

		static int ENUM_RENDER_NONE;
		static int ENUM_RENDER_VELOCITY;
		static int ENUM_RENDER_ANGULAR_VELOCITY;
		static int ENUM_RENDER_DENSITY;

		SimulatorBase();
		virtual ~SimulatorBase();

		void init(int argc, char **argv, const char *simName);
		void buildModel();
		void cleanup();

		void renderFluid(const unsigned int fluidModelIndex, float *fluidColor);

		void readParameters();
		void partioExport();
		void writeParticles(const std::string &fileName, FluidModel *model);
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
		void pointShaderBegin(Shader *shader, const float *col, const Real minVal, const Real maxVal, const bool useTexture = false, float const* color_map = nullptr);
		void pointShaderEnd(Shader *shader, const bool useTexture = false);
		Utilities::SceneLoader::Scene& getScene() { return m_scene; }

		std::vector<std::vector<unsigned int>>& getSelectedParticles() { return m_selectedParticles; }
		bool getUseParticleCaching() const { return m_useParticleCaching; }
		void setUseParticleCaching(bool val) { m_useParticleCaching = val; }

		const std::string& getColorField(const unsigned int fluidModelIndex) {	return m_colorField[fluidModelIndex]; }
		void setColorField(const unsigned int fluidModelIndex, const std::string& fieldName) { m_colorField[fluidModelIndex] = fieldName; }

		int getColorMapType(const unsigned int fluidModelIndex) const { return m_colorMapType[fluidModelIndex]; }
		void setColorMapType(const unsigned int fluidModelIndex, int val) { m_colorMapType[fluidModelIndex] = val; }
		Real getRenderMaxValue(const unsigned int fluidModelIndex) const { return m_renderMaxValue[fluidModelIndex]; }
		void setRenderMaxValue(const unsigned int fluidModelIndex, Real val) { m_renderMaxValue[fluidModelIndex] = val; }
		Real getRenderMinValue(const unsigned int fluidModelIndex) const { return m_renderMinValue[fluidModelIndex]; }
		void setRenderMinValue(const unsigned int fluidModelIndex, Real val) { m_renderMinValue[fluidModelIndex] = val; }
	};
}
 
#endif