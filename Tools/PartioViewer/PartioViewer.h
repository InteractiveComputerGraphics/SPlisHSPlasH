#pragma once

#include <string>
#include "SPlisHSPlasH/TriangleMesh.h"
#include "extern/partio/src/lib/Partio.h"
#include "GUI/PartioViewer_GUI_Base.h"
#include "extern/AntTweakBar/include/AntTweakBar.h"

namespace SPH
{
	struct Fluid
	{
		Partio::ParticlesDataMutable* partioData;
		std::vector<unsigned int> selectedParticles;
		std::vector<unsigned int> visibleParticles;
		std::string inputFile;
		std::string currentFile;
		unsigned int posIndex;
		unsigned int m_colorField;
		float m_renderMinValue;
		float m_renderMaxValue;
		unsigned int m_colorMapType;
	};

	struct Boundary
	{
		TriangleMesh mesh;
		std::vector<Vector3r> x0;
		Vector3f t;
		Matrix3f R;
		bool isWall;
		Vector4f color;
	};

	class PartioViewer
	{
	protected:
		int m_frameIndex;
		int m_startFrame;
		int m_endFrame;
		int m_firstRBIndex;
		unsigned int m_numberOfRigidBodies;
		Real m_particleRadius;
		Eigen::Vector3f m_planeNormal;
		Eigen::Vector3f m_planePoint;
		bool m_renderSequence;
		bool m_renderVideo;
		bool m_overWrite;
		std::string m_ffmpegPath;
		std::string m_rbDataFile;
		FILE *m_ffmpegPipe;
		unsigned int m_width;
		unsigned int m_height;
		unsigned int m_fps;
		static FILE *m_jpegFile;
		bool m_doPause;
		std::vector<Fluid> m_fluids;
		std::vector<Boundary> m_boundaries;		
		std::string m_exePath, m_outPath;		
		bool m_usePlane;
		bool m_useRBData;
		Eigen::AlignedBox3f m_fluidBoundingBox;	
		std::vector<unsigned char> m_image;
		PartioViewer_GUI_Base *m_gui;
		int m_argc;
		char **m_argv;
		std::string m_defaultColorFieldName;
		float m_defaultRenderMinValue;
		float m_defaultRenderMaxValue;
		unsigned int m_defaultColorMapType;

		void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale);
		bool readPartioFile(const std::string &fileName, Partio::ParticlesDataMutable* &partioData, unsigned int &posIndex, const bool printInfo = false);
		bool readRigidBodyData(std::string fileName, const bool first = false);
		unsigned int getFrameIndexFromFile(const std::string &inputFileName);
		std::string convertFileName(const std::string &inputFileName, const std::string currentFrame);
		std::string zeroPadding(const unsigned int number, const unsigned int length);
		bool imagesExist(const unsigned int frameIndex);
		void addImageDataToVideo();
		void saveImageData();
		static void outputJpg(unsigned char oneByte);
		static void saveImage(const std::string &fileName, const void *pixels, unsigned int width, unsigned int height);
		void initGUI();

	public:
		PartioViewer();
		~PartioViewer();

		int run(int argc, char **argv);
		bool setFrame(const unsigned int index);
		void reset();
		void updateBoundingBox();
		bool updateData();

		void particleInfo();
		bool nextFrame();
		bool prevFrame();
		void saveFrame();
		void generateSequence();
		void generateVideo();
		void timeStep();
		
		int getFrameIndex() const { return m_frameIndex; }
		bool getUsePlane() const { return m_usePlane; }
		void setUsePlane(bool val);
		Eigen::Vector3f getPlaneNormal() const { return m_planeNormal; }
		void setPlaneNormal(const Eigen::Vector3f &val);
		Eigen::Vector3f getPlanePoint() const { return m_planePoint; }
		void setPlanePoint(const Eigen::Vector3f &val);
		std::vector<SPH::Fluid>& getFluids() { return m_fluids; }
		std::vector<SPH::Boundary>& getBoundaries() { return m_boundaries; }
		Real getParticleRadius() const { return m_particleRadius; }
		void setParticleRadius(Real val) { m_particleRadius = val; }
		Eigen::AlignedBox3f getFluidBoundingBox() const { return m_fluidBoundingBox; }
		void setFluidBoundingBox(const Eigen::AlignedBox3f &val) { m_fluidBoundingBox = val; }
		std::string getExePath() const { return m_exePath; }
		void setExePath(const std::string &val) { m_exePath = val; }
		void setStartFrame(const int val) { m_startFrame = val; }
		int getStartFrame() const { return m_startFrame; }
		void setEndFrame(const int val) { m_endFrame = val; }
		int getEndFrame() const { return m_endFrame; }
		void setFPS(const int val) { m_fps = val; }
		int getFPS() const { return m_fps; }
		int getArgc() const { return m_argc; }
		char ** getArgv() const { return m_argv; }
		unsigned int getWidth() const { return m_width; }
		unsigned int getHeight() const { return m_height; }
		bool getPause() const { return m_doPause; }
		void setPause(bool val) { m_doPause = val; }
	};
}
