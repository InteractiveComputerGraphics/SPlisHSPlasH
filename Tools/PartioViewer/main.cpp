#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include <iostream>
#include "GL/glew.h"
#include "Visualization/MiniGL.h"
#include "GL/glut.h"
#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/PoissonDiskSampling.h"
#include "Utilities/FileSystem.h"
#include <cfloat>
#include "Utilities/Version.h"
#include "Visualization/Selection.h"
#include "extern/partio/src/lib/Partio.h"
#include "Visualization/colormaps/colormap_jet.h"
#include "Visualization/colormaps/colormap_plasma.h"
#include "extern/cxxopts/cxxopts.hpp"
#include "GL/freeglut_ext.h"
#include <regex>
#include <stdio.h>
#include "extern/toojpeg/toojpeg.h"
#include "Utilities/BinaryFileReaderWriter.h"


// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

INIT_TIMING
INIT_LOGGING

using namespace SPH;
using namespace std;
using namespace Utilities;

struct Fluid
{
	Partio::ParticlesDataMutable* partioData;
	std::vector<unsigned int> selectedParticles;
	std::vector<unsigned int> visibleParticles;
	string inputFile;
	string currentFile;
	unsigned int posIndex;
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

void initShaders();
void render();
void render(const Fluid &fluid, float *fluidColor);
void render(const Boundary &boundary, float *boundaryColor);
void pointShaderBegin(Shader *shader, const float *col, const Real minVal, const Real maxVal, const bool useTexture = false, float const* color_map = nullptr);
void pointShaderEnd(Shader *shader, const bool useTexture = false);
void timeStep();
void updateBoundingBox();
void selection(const Vector2i &start, const Vector2i &end, void *clientData);
void particleInfo();
void saveFrame();
bool nextFrame();
bool prevFrame();
bool setFrame(const unsigned int index);
bool updateData();
void saveImage(const std::string &fileName);
void saveImageData();
void addImageDataToVideo();
void getImage();
void initGUI();
void generateVideo();
void generateSequence();
void flipImage(int width, int height, unsigned char *image);
bool readPartioFile(const std::string &fileName, Partio::ParticlesDataMutable* &partioData, unsigned int &posIndex, const bool printInfo = false);
bool readRigidBodyData(std::string fileName, const bool first = false);
bool imagesExist(const unsigned int frameIndex);
unsigned int getFrameIndexFromFile(const std::string &inputFileName);
std::string convertFileName(const std::string &inputFileName, const std::string currentFrame);
void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale);
void TW_CALL setFrameIndex(const void *value, void *clientData);
void TW_CALL getFrameIndex(void *value, void *clientData);
void TW_CALL setShowBBox(const void *value, void *clientData);
void TW_CALL getShowBBox(void *value, void *clientData);
void TW_CALL setUsePlane(const void *value, void *clientData);
void TW_CALL getUsePlane(void *value, void *clientData);
void TW_CALL setPlaneNormal(const void *value, void *clientData);
void TW_CALL getPlaneNormal(void *value, void *clientData);
void TW_CALL setPlanePoint(const void *value, void *clientData);
void TW_CALL getPlanePoint(void *value, void *clientData);

FILE *jpegFile;
bool doPause = true;
std::vector<Fluid> fluids;
std::vector<Boundary> boundaries;
unsigned int numberOfRigidBodies = 0;
int firstRBIndex = 1;
string exePath, dataPath, outPath;
Real particleRadius = 0.025;
Vector3f planeNormal(0, 0, -1);
Vector3f planePoint(0, 0, 0);
bool usePlane = false;
bool useRBData = false;
Eigen::AlignedBox3f fluidBoundingBox;
Shader shader_vector;
Shader shader_scalar;
Shader shader_vector_map;
Shader shader_scalar_map;
Shader meshShader;
GLuint textureMap;
GLint context_major_version;
GLint context_minor_version;
unsigned int colorField = 0;
std::string colorFieldName = "velocity";
unsigned int colorMapType = 1;
float renderMinValue = 0.0;
float renderMaxValue = 10.0;
int frameIndex = 1;
int startFrame = -1;
int endFrame = -1;
bool renderSequence = false;
bool renderVideo = false;
bool overWrite = true;
std::string ffmpegPath = "ffmpeg";
std::string rbDataFile = "";
bool renderWalls = false;
bool showBBox = false;
Vector3r camPos(0.0, 3.0, 10.0);
Vector3r camLookat(0, 0, 0);
FILE *ffmpegPipe = nullptr;
unsigned int width = 1280;
unsigned int height = 960;
unsigned int fps = 25;
std::vector<unsigned char> image;


std::istream& operator >> (std::istream& istream, Vector3r& v)
{
	return istream >> std::skipws >> v[0] >> v[1] >> v[2];
}


// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	Utilities::logger.addSink(unique_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));

	LOG_INFO << "Git refspec: " << GIT_REFSPEC;
	LOG_INFO << "Git SHA1: " << GIT_SHA1;
	LOG_INFO << "Git status: " << GIT_LOCAL_STATUS;

	exePath = FileSystem::getProgramPath();
	dataPath = FileSystem::normalizePath(exePath + "/" + std::string(SPH_DATA_PATH));
	try
	{
		cxxopts::Options options(argv[0], "PartioViewer - Visualize partio data.");
		options
			.positional_help("[partio file]")
			.show_positional_help();

		options.add_options()
			("h,help", "Print help")
			("renderSequence", "Render a sequence from startFrame to endFrame as jpeg.")
			("renderVideo", "Render a sequence from startFrame to endFrame as video."
							"This function requires ffmpeg which must be in the PATH or the ffmpegPath parameter must be set.")
			("noOverwrite", "Do not overwrite existing frames when using --renderSequence option."
							"Existing frames are not loaded at all which accelerates the image sequence generation.")
			("o,outdir", "Output directory for images", cxxopts::value<std::string>())
			("rbData", "Rigid body data to visualize (bin file)", cxxopts::value<std::string>())
			("ffmpegPath", "Path of the ffmpeg excutable.", cxxopts::value<std::string>())
			("width", "Width of the image in pixels.", cxxopts::value<int>()->default_value("1024"))
			("height", "Height of the image in pixels.", cxxopts::value<int>()->default_value("768"))
			("fps", "Frame rate of video.", cxxopts::value<int>()->default_value("25"))
			("r,radius", "Particle radius", cxxopts::value<Real>()->default_value("0.025"))
			("s,startFrame", "Start frame (only used if value is >= 0)", cxxopts::value<int>()->default_value("-1"))
			("e,endFrame", "End frame (only used if value is >= 0)", cxxopts::value<int>()->default_value("-1"))
			("colorField", "Name of field that is used for the color.", cxxopts::value<std::string>()->default_value("velocity"))
			("colorMapType", "Color map (0=None, 1=Jet, 2=Plasma)", cxxopts::value<unsigned int>()->default_value("1"))
			("renderMinValue", "Min value of field.", cxxopts::value<float>()->default_value("0.0"))
			("renderMaxValue", "Max value of field.", cxxopts::value<float>()->default_value("10.0"))
			("camPos", "Camera position (e.g. --camPos \"0 1 5\")", cxxopts::value<Vector3r>()->default_value("0 3 10"))
			("camLookat", "Camera lookat (e.g. --camLookat \"0 0 0\")", cxxopts::value<Vector3r>()->default_value("0 0 0"))
			;

		options.add_options("invisible")
			("partio-file", "Partio file", cxxopts::value<std::vector<std::string>>());

		options.parse_positional({ "partio-file" });
		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			LOG_INFO << options.help({ "", "Group" });
			exit(0);
		}

		if (result.count("partio-file"))
		{
			std::vector<std::string> inputFile = result["partio-file"].as<std::vector<std::string>>();

			// if a directory is given as input, assume that the partio files are in directory partio
			// and rigid body data is in rigid_bodies
			if (inputFile.size() == 1)
			{
				// check if argument is a directory
				if (FileSystem::isDirectory(inputFile[0]))
				{
					std::string oDir = inputFile[0];
					std::string partioDir = FileSystem::normalizePath(inputFile[0] + "/partio");

					// check if partio directory exists
					if (FileSystem::isDirectory(partioDir))
					{
						// get all files in directory
						std::vector<std::string> fileNames;
						FileSystem::getFilesInDirectory(partioDir, fileNames);
						std::cout << "test\n";
						std::vector<std::string> fileBases;
					
						// filter *.bgeo files
						fileNames.erase(std::remove_if(fileNames.begin(), fileNames.end(), [](const std::string &s) { return FileSystem::getFileExt(s) != "bgeo"; }), fileNames.end());

						if (fileNames.size() > 0)
						{
							// create list of pairs <index, filename>
							std::vector<std::pair<unsigned int, std::string>> framePairs;
							framePairs.resize(fileNames.size());
							for (int i = 0; i < fileNames.size(); i++)
								framePairs[i] = { getFrameIndexFromFile(fileNames[i]), fileNames[i] };

							// sort list with respect to index
							std::sort(framePairs.begin(), framePairs.end());

							while (framePairs.size() > 0)
							{
								// found a new fluid phase
								fileBases.push_back(framePairs[0].second);

								// remove all other file names of this phase from the list
								std::string test = convertFileName(framePairs[0].second, "#");
								framePairs.erase(std::remove_if(framePairs.begin(), framePairs.end(), [&](const std::pair<unsigned int, std::string> &s) { return test == convertFileName(s.second, "#"); }), framePairs.end());
							}

							// generate list of input files
							inputFile.clear();
							inputFile.resize(fileBases.size());
							for (int i = 0; i < fileBases.size(); i++)
							{
								inputFile[i] = FileSystem::normalizePath(partioDir + "/" + fileBases[i]);
							}
						}
					}

					// only search for files if user did not define the file 
					if (!result.count("rbData"))
					{
						// check if directory exists
						std::string rbDir = FileSystem::normalizePath(oDir + "/rigid_bodies");
						if (FileSystem::isDirectory(rbDir))
						{
							// get all files in directory
							std::vector<std::string> fileNames;
							FileSystem::getFilesInDirectory(rbDir, fileNames);

							// filter *.bin files
							fileNames.erase(std::remove_if(fileNames.begin(), fileNames.end(), [](const std::string &s) { return FileSystem::getFileExt(s) != "bin"; }), fileNames.end());

							if (fileNames.size() > 0)
							{
								// create list of pairs <index, filename>
								std::vector<std::pair<unsigned int, std::string>> framePairs;
								framePairs.resize(fileNames.size());
								for (int i = 0; i < fileNames.size(); i++)
									framePairs[i] = { getFrameIndexFromFile(fileNames[i]), fileNames[i] };

								// sort list with respect to index
								std::sort(framePairs.begin(), framePairs.end());

								// get file with lowest index
								rbDataFile = FileSystem::normalizePath(rbDir + "/" + framePairs[0].second);
								useRBData = true;
							}
						}
					}
				}
			}


			fluids.resize(inputFile.size());
			for (int i=0; i < inputFile.size(); i++)
			{
				fluids[i].inputFile = inputFile[i];
				if (FileSystem::isRelativePath(inputFile[i]))
					fluids[i].inputFile = FileSystem::normalizePath(exePath + "/" + inputFile[i]);
				fluids[i].currentFile = fluids[i].inputFile;
			}
		}
		else
		{
			LOG_INFO << options.help({ "", "Group" });
			exit(0);
		}

		outPath = exePath;
		if (result.count("outdir"))
			outPath = result["outdir"].as<std::string>();

		if (result.count("rbData"))
		{
			rbDataFile = result["rbData"].as<std::string>();
			useRBData = true;
		}

		if (result.count("ffmpegPath"))
			ffmpegPath = result["ffmpegPath"].as<std::string>();

		particleRadius = 0.025;
		if (result.count("radius"))
			particleRadius = result["radius"].as<Real>();
		LOG_INFO << "Radius: " << particleRadius;

		if (result.count("startFrame"))
			startFrame = result["startFrame"].as<int>();
		LOG_INFO << "Start frame: " << startFrame;

		if (result.count("endFrame"))
			endFrame = result["endFrame"].as<int>();
		LOG_INFO << "End frame: " << endFrame;

		if (result.count("renderSequence"))
		{
			renderSequence = true;
			if (result.count("noOverwrite"))
				overWrite = false;
		}

		if (result.count("renderVideo"))
			renderVideo = true;

		renderMinValue = 0.0;
		if (result.count("renderMinValue"))
			renderMinValue = result["renderMinValue"].as<float>();

		renderMaxValue = 10.0;
		if (result.count("renderMaxValue"))
			renderMaxValue = result["renderMaxValue"].as<float>();

		colorMapType = 1;
		if (result.count("colorMapType"))
			colorMapType = result["colorMapType"].as<unsigned int>();

		colorFieldName = "velocity";
		if (result.count("colorField"))
			colorFieldName = result["colorField"].as<std::string>();

		if (result.count("camPos"))
			camPos = result["camPos"].as<Vector3r>();

		if (result.count("camLookat"))
			camLookat = result["camLookat"].as<Vector3r>();

		if (result.count("width"))
			width = result["width"].as<int>();

		if (result.count("height"))
			height = result["height"].as<int>();

		if (result.count("fps"))
			fps = result["fps"].as<int>();

		image.resize(3 * width * height);
	}
	catch (const cxxopts::OptionException& e)
	{
		LOG_INFO << "error parsing options: " << e.what();
		exit(1);
	}

	frameIndex = getFrameIndexFromFile(fluids[0].inputFile);
	if (frameIndex == 0xffffffff)
		frameIndex = 0;
	firstRBIndex = frameIndex;

	// read particle data
	if (!updateData())
	{
		LOG_ERR << "Unable to read partio file or did not find position data.";
			exit(1);
	}

	// OpenGL
	MiniGL::init(argc, argv, width, height, 0, 0, "Partio Viewer");
	MiniGL::initLights();
	MiniGL::initTweakBarParameters();
	MiniGL::getOpenGLVersion(context_major_version, context_minor_version);
	MiniGL::setViewport(40.0, 0.1f, 500.0, camPos, camLookat);
	MiniGL::setSelectionFunc(selection, nullptr);
	MiniGL::addKeyFunc('i', particleInfo);
	MiniGL::addKeyFunc('+', nextFrame);
	MiniGL::addKeyFunc('-', prevFrame);
	MiniGL::addKeyFunc('s', saveFrame);
	MiniGL::addKeyFunc('v', generateVideo);
	MiniGL::addKeyFunc('j', generateSequence);
	MiniGL::addKeyFunc(' ', [&] { doPause = !doPause; });
	MiniGL::addKeyFunc('r', [&] { setFrameIndex(&firstRBIndex, nullptr); TwRefreshBar(MiniGL::getTweakBar());  });

	initGUI();
	if (showBBox)
		updateBoundingBox();

	if (MiniGL::checkOpenGLVersion(3, 3))
		initShaders();
	
	MiniGL::setClientSceneFunc(render);
	MiniGL::setClientIdleFunc(1, timeStep);

	if ((renderSequence || renderVideo) && (startFrame >= 0))
		setFrame(startFrame);

	glutMainLoop();

	for (size_t i = 0; i < fluids.size(); i++)
		fluids[i].partioData->release();
	Timing::printAverageTimes();
	Timing::printTimeSums();
	
	return 0;
}

void timeStep()
{
	if (!renderSequence && !renderVideo)
	{
		if (!doPause)
		{
			if (!nextFrame())
				doPause = true;
		}
	}
	else if (renderVideo)
	{
		generateVideo();
		glutLeaveMainLoop();
	}
	else if (renderSequence)
	{
		generateSequence();
		glutLeaveMainLoop();
	}
}

void generateSequence()
{
	width = MiniGL::getWidth();
	height = MiniGL::getHeight();

	// height must be divisible by 2
	if (height % 2 == 1)
		height--;
	// width must be divisible by 2
	if (width % 2 == 1)
		width--;

	image.resize(3 * width*height);

	if ((startFrame >= 0) && (!imagesExist(startFrame)))
		setFrame(startFrame);

	bool chk = true;
	while (chk && ((frameIndex <= endFrame) || (endFrame < 0)))
	{
		if (overWrite || !imagesExist(frameIndex))
		{
			if (!updateData())
			{
				frameIndex--;
				chk = false;
				break;
			}

			MiniGL::viewport();
			render();
			glutSwapBuffers();

			getImage();
			saveImageData();
		}
		frameIndex++;
	}
}


void generateVideo()
{
	width = MiniGL::getWidth();
	height = MiniGL::getHeight();

	// height must be divisible by 2
	if (height % 2 == 1)
		height--;
	// width must be divisible by 2
	if (width % 2 == 1)
		width--;

	// open pipe to ffmpeg
	stringstream sstm;

	FileSystem::makeDirs(outPath);
	std::string videoFile = FileSystem::normalizePath(outPath + "/output.mp4");
	sstm << ffmpegPath << " -y -hide_banner -nostats -loglevel panic -r " << fps << " -f rawvideo -vcodec rawvideo -pix_fmt rgb24 -s "
		<< width << "x" << height << " -i - -threads 0 -preset fast -pix_fmt yuv420p -crf 21 " << videoFile << " > ffmpeg.log";

#ifdef WIN32
	if (!(ffmpegPipe = _popen(sstm.str().c_str(), "wb")))
#else
	if (!(ffmpegPipe = popen(sstm.str().c_str(), "w")))
#endif
	{
		LOG_ERR << "Cannot open pipe to ffmpeg.";
		return;
	}

	image.resize(3 * width*height);

	if (startFrame >= 0)
		setFrame(startFrame);

	bool chk = true;
	while (chk && ((frameIndex <= endFrame) || (endFrame < 0)))
	{
		MiniGL::viewport();
		render();
		glutSwapBuffers();

		getImage();
		addImageDataToVideo();
		chk = nextFrame();
	}
	fflush(ffmpegPipe);
	fclose(ffmpegPipe);
	ffmpegPipe = nullptr;
	LOG_INFO << "Generated video: " << videoFile;
}

void initGUI()
{
	TwAddVarCB(MiniGL::getTweakBar(), "frameIndex", TW_TYPE_UINT32, setFrameIndex, getFrameIndex, nullptr, " label='Frame index' min=0 group=General");
	TwAddVarRW(MiniGL::getTweakBar(), "particleRadius", TW_TYPE_REAL, &particleRadius, " label='Particle radius' min=0.001 group=General");

	TwAddVarRW(MiniGL::getTweakBar(), "renderWalls", TW_TYPE_BOOL32, &renderWalls, " label='Render walls' group=Visualization");
	TwAddVarCB(MiniGL::getTweakBar(), "showBBox", TW_TYPE_BOOL32, setShowBBox, getShowBBox, nullptr, " label='Show bounding box' group=Visualization");

	TwAddVarCB(MiniGL::getTweakBar(), "usePlane", TW_TYPE_BOOL32, setUsePlane, getUsePlane, nullptr, " label='Use plane' group=Visualization");
	TwAddVarCB(MiniGL::getTweakBar(), "planePoint", TW_TYPE_DIR3F, setPlanePoint, getPlanePoint, nullptr, " label='Plane point' group=Visualization");
	TwAddVarCB(MiniGL::getTweakBar(), "planeNormal", TW_TYPE_DIR3F, setPlaneNormal, getPlaneNormal, nullptr, " label='Plane normal' group=Visualization");
	

	TwType enumType = TwDefineEnum("ColorFieldEnum", NULL, 0);
	std::ostringstream oss;
	int idx = 0;
	int velIdx = 0;

	if (fluids[0].partioData)
	{
		bool found = false;
		for (int i = 0; i < fluids[0].partioData->numAttributes(); i++)
		{
			Partio::ParticleAttribute attr;
			fluids[0].partioData->attributeInfo(i, attr);
			if ((attr.type == Partio::FLOAT) || (attr.type == Partio::VECTOR))
			{
				if (idx != 0)
					oss << ", ";
				oss << idx << " {" << attr.name.c_str() << "}";
				
				if (attr.name == colorFieldName)
				{
					colorField = idx;
					found = true;
				}
				if (attr.name == "velocity")
					velIdx = idx;
			}
			idx++;
		}
		// Choose velocity as default
		if (!found)
			colorField = velIdx;
	}
	
	std::string enumStr = " label='Color field' enum='" + oss.str() + "' group='Visualization'";
	enumStr = enumStr + " help='Choose vector or scalar field for particle coloring.'";
	TwAddVarRW(MiniGL::getTweakBar(), "ColorField", enumType, &colorField, enumStr.c_str());

	TwType enumType2 = TwDefineEnum("ColorMapTypeEnum", NULL, 0);
	std::string str = " label='Color map' enum='0 {None}, 1 {Jet}, 2 {Plasma}' group='Visualization' help='Choose a color map.'";
	TwAddVarRW(MiniGL::getTweakBar(), "ColorMapType", enumType2, &colorMapType, str.c_str());
	str = " label='Min. value (shader)' step=0.001 precision=3 group='Visualization' help='Minimal value used for color-coding the color field in the rendering process.'";
	TwAddVarRW(MiniGL::getTweakBar(), "RenderMinValue", TW_TYPE_FLOAT, &renderMinValue, str.c_str());
	str = " label='Max. value (shader)' step=0.001 precision=3 group='Visualization' help='Maximal value used for color-coding the color field in the rendering process.'";
	TwAddVarRW(MiniGL::getTweakBar(), "RenderMaxValue", TW_TYPE_FLOAT, &renderMaxValue, str.c_str());

	TwAddVarRW(MiniGL::getTweakBar(), "StartFrame", TW_TYPE_INT32, &startFrame, " label='Start frame' min=0 group=Export");
	TwAddVarRW(MiniGL::getTweakBar(), "EndFrame", TW_TYPE_INT32, &endFrame, " label='End frame' min=0 group=Export");
	TwAddVarRW(MiniGL::getTweakBar(), "FPS", TW_TYPE_INT32, &fps, " label='FPS' min=1 group=Export");
}

void initShaders()
{
	string vertFile = exePath + "/resources/shaders/vs_points_vector.glsl";
	string fragFile = exePath + "/resources/shaders/fs_points.glsl";
	shader_vector.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	shader_vector.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	shader_vector.createAndLinkProgram();
	shader_vector.begin();
	shader_vector.addUniform("modelview_matrix");
	shader_vector.addUniform("projection_matrix");
	shader_vector.addUniform("radius");
	shader_vector.addUniform("viewport_width");
	shader_vector.addUniform("color");
	shader_vector.addUniform("min_scalar");
	shader_vector.addUniform("max_scalar");
	shader_vector.end();

	string vertFileScalar = exePath + "/resources/shaders/vs_points_scalar.glsl";
	shader_scalar.compileShaderFile(GL_VERTEX_SHADER, vertFileScalar);
	shader_scalar.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	shader_scalar.createAndLinkProgram();
	shader_scalar.begin();
	shader_scalar.addUniform("modelview_matrix");
	shader_scalar.addUniform("projection_matrix");
	shader_scalar.addUniform("radius");
	shader_scalar.addUniform("viewport_width");
	shader_scalar.addUniform("color");
	shader_scalar.addUniform("min_scalar");
	shader_scalar.addUniform("max_scalar");
	shader_scalar.end();

	string fragFileMap = exePath + "/resources/shaders/fs_points_colormap.glsl";
	shader_vector_map.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	shader_vector_map.compileShaderFile(GL_FRAGMENT_SHADER, fragFileMap);
	shader_vector_map.createAndLinkProgram();
	shader_vector_map.begin();
	shader_vector_map.addUniform("modelview_matrix");
	shader_vector_map.addUniform("projection_matrix");
	shader_vector_map.addUniform("radius");
	shader_vector_map.addUniform("viewport_width");
	shader_vector_map.addUniform("color");
	shader_vector_map.addUniform("min_scalar");
	shader_vector_map.addUniform("max_scalar");
	shader_vector_map.end();

	shader_scalar_map.compileShaderFile(GL_VERTEX_SHADER, vertFileScalar);
	shader_scalar_map.compileShaderFile(GL_FRAGMENT_SHADER, fragFileMap);
	shader_scalar_map.createAndLinkProgram();
	shader_scalar_map.begin();
	shader_scalar_map.addUniform("modelview_matrix");
	shader_scalar_map.addUniform("projection_matrix");
	shader_scalar_map.addUniform("radius");
	shader_scalar_map.addUniform("viewport_width");
	shader_scalar_map.addUniform("color");
	shader_scalar_map.addUniform("min_scalar");
	shader_scalar_map.addUniform("max_scalar");
	shader_scalar_map.end();

	vertFile = exePath + "/resources/shaders/vs_smooth.glsl";
	fragFile = exePath + "/resources/shaders/fs_smooth.glsl";
	meshShader.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	meshShader.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	meshShader.createAndLinkProgram();
	meshShader.begin();
	meshShader.addUniform("modelview_matrix");
	meshShader.addUniform("projection_matrix");
	meshShader.addUniform("surface_color");
	meshShader.addUniform("shininess");
	meshShader.addUniform("specular_factor");
	meshShader.end();

	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &textureMap);
}

void renderAABB(const Eigen::AlignedBox3f &aabb, float *color)
{
	const Vector3r a = aabb.corner(Eigen::AlignedBox3f::BottomLeftFloor).cast<Real>();
	const Vector3r b = aabb.corner(Eigen::AlignedBox3f::BottomRightFloor).cast<Real>();
	const Vector3r c = aabb.corner(Eigen::AlignedBox3f::TopRightFloor).cast<Real>();
	const Vector3r d = aabb.corner(Eigen::AlignedBox3f::TopLeftFloor).cast<Real>();
	const Vector3r e = aabb.corner(Eigen::AlignedBox3f::BottomLeftCeil).cast<Real>();
	const Vector3r f = aabb.corner(Eigen::AlignedBox3f::BottomRightCeil).cast<Real>();
	const Vector3r g = aabb.corner(Eigen::AlignedBox3f::TopRightCeil).cast<Real>();
	const Vector3r h = aabb.corner(Eigen::AlignedBox3f::TopLeftCeil).cast<Real>();

	const float w = 1.0;
	MiniGL::drawVector(a, b, w, color);
	MiniGL::drawVector(b, c, w, color);
	MiniGL::drawVector(c, d, w, color);
	MiniGL::drawVector(d, a, w, color);

	MiniGL::drawVector(e, f, w, color);
	MiniGL::drawVector(f, g, w, color);
	MiniGL::drawVector(g, h, w, color);
	MiniGL::drawVector(h, e, w, color);

	MiniGL::drawVector(a, e, w, color);
	MiniGL::drawVector(b, f, w, color);
	MiniGL::drawVector(c, g, w, color);
	MiniGL::drawVector(d, h, w, color);
}

void render()
{
	for (size_t i = 0; i < fluids.size(); i++)
	{
		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		MiniGL::hsvToRgb(0.61f - 0.1f*i, 0.66f, 0.9f, fluidColor);
		render(fluids[i], fluidColor);
	}

	float boundaryColor[4] = { 0.4f, 0.4f, 0.4f, 1.0f };
	for (size_t i = 0; i < boundaries.size(); i++)
	{
		render(boundaries[i], boundaryColor);
	}
}

void render(const Fluid &fluid, float *fluidColor)
{
	MiniGL::coordinateSystem();

	float gridColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
	MiniGL::drawGrid_xz(gridColor);

	// Draw simulation model
	const unsigned int nParticles = (unsigned int)fluid.partioData->numParticles();

	Partio::ParticleAttribute posAttr;
	fluid.partioData->attributeInfo(fluid.posIndex, posAttr);
	const float* partioX = fluid.partioData->data<float>(posAttr, 0);

	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		Shader *shader_vec = &shader_vector_map;
		Shader *shader_s = &shader_scalar_map;
		float const *color_map = nullptr;
		if (colorMapType == 1)
			color_map = reinterpret_cast<float const*>(colormap_jet);
		else if (colorMapType == 2)
			color_map = reinterpret_cast<float const*>(colormap_plasma);

		if (colorMapType == 0)
		{
			shader_vec = &shader_vector;
			shader_s = &shader_scalar;
		}

		if (fluid.partioData->numAttributes() == 0)
			pointShaderBegin(shader_s, fluidColor, renderMinValue, renderMaxValue, false);
		else 
		{
			Partio::ParticleAttribute attr;
			fluid.partioData->attributeInfo(colorField, attr);

			if (attr.type == Partio::VECTOR)
				pointShaderBegin(shader_vec, fluidColor, renderMinValue, renderMaxValue, true, color_map);
			else if (attr.type == Partio::FLOAT)
				pointShaderBegin(shader_s, fluidColor, renderMinValue, renderMaxValue, true, color_map);
		}


		if (nParticles > 0)
		{
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, partioX);

			if (fluid.partioData->numAttributes() > 0)
			{
				Partio::ParticleAttribute attr;
				fluid.partioData->attributeInfo(colorField, attr);

				const float* partioVals = NULL;
				if (attr.type == Partio::VECTOR)
				{
					glEnableVertexAttribArray(1);
					partioVals = fluid.partioData->data<float>(attr, 0);
					glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, &partioVals[0]);
				}
				else if (attr.type == Partio::FLOAT)
				{
					glEnableVertexAttribArray(1);
					partioVals = fluid.partioData->data<float>(attr, 0);
					glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, &partioVals[0]);
				}
			}

			if (usePlane)
				glDrawElements(GL_POINTS, (GLsizei)fluid.visibleParticles.size(), GL_UNSIGNED_INT, fluid.visibleParticles.data());
			else
				glDrawArrays(GL_POINTS, 0, nParticles);

			glDisableVertexAttribArray(0);
			//			glDisableVertexAttribArray(1);
		}

		if (fluid.partioData->numAttributes() == 0)
			pointShaderEnd(shader_s, false);
		else 
		{
			Partio::ParticleAttribute attr;
			fluid.partioData->attributeInfo(colorField, attr);
			if (attr.type == Partio::VECTOR)
				pointShaderEnd(shader_vec, true);
			else if (attr.type == Partio::FLOAT)
				pointShaderEnd(shader_s, true);
		}
	}
	else
	{
		const Real supportRadius = particleRadius*4.0;
		float fluidColor[4] = { 0.1f, 0.2f, 0.6f, 1.0f };

		Partio::ParticleAttribute attr;
		if (fluid.partioData->numAttributes() > 0)
			fluid.partioData->attributeInfo(colorField, attr);

		glPointSize(4.0);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < nParticles; i++)
		{
			if (fluid.partioData->numAttributes() > 0)
			{
				float v = 0.0;
				const float* partioVals = NULL;
				if (attr.type == Partio::VECTOR)
				{
					partioVals = fluid.partioData->data<float>(attr, 0);
					v = sqrt(partioVals[3 * i] * partioVals[3 * i] + partioVals[3 * i + 1] * partioVals[3 * i + 1] + partioVals[3 * i + 2] * partioVals[3 * i + 2]);
				}
				else if (attr.type == Partio::FLOAT)
				{
					partioVals = fluid.partioData->data<float>(attr, 0);
					v = partioVals[3 * i];
				}

				v = 0.5f*((v - renderMinValue) / (renderMaxValue - renderMinValue));
				v = min(128.0f*v*v, 0.5f);
				float fluidColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
				MiniGL::hsvToRgb(0.55f, 1.0f, 0.5f + v, fluidColor);
				glColor3fv(fluidColor);
			}
			else 
				glColor3fv(fluidColor);
			glVertex3fv(&partioX[3*i]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}

	float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		pointShaderBegin(&shader_vector, &red[0], renderMinValue, renderMaxValue);
		if (fluid.selectedParticles.size() > 0)
		{
			glUniform1f(shader_vector.getUniform("radius"), (float)particleRadius*1.05f);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FALSE, GL_FALSE, 0, partioX);

			if (fluid.partioData->numAttributes() > 0)
			{
				Partio::ParticleAttribute attr;
				fluid.partioData->attributeInfo(colorField, attr);
				const float* partioVals = NULL;
				if (attr.type == Partio::VECTOR)
				{
					glEnableVertexAttribArray(1);
					partioVals = fluid.partioData->data<float>(attr, 0);
					glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, &partioVals[0]);
				}
				else if (attr.type == Partio::FLOAT)
				{
					glEnableVertexAttribArray(1);
					partioVals = fluid.partioData->data<float>(attr, 0);
					glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, &partioVals[0]);
				}
			}

			glDrawElements(GL_POINTS, (GLsizei)fluid.selectedParticles.size(), GL_UNSIGNED_INT, fluid.selectedParticles.data());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}
		pointShaderEnd(&shader_vector);
	}
	else
	{
		if (fluid.selectedParticles.size() > 0)
		{
			glPointSize(4.0);
			glDisable(GL_LIGHTING);
			glBegin(GL_POINTS);
			for (unsigned int i = 0; i < fluid.selectedParticles.size(); i++)
			{
				glColor3fv(red);
				glVertex3fv(&partioX[3* fluid.selectedParticles[i]]);
			}
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}

	// Render bounding box - fluid
	float col[4] = { 0.8,0.8,0.8,1 };
	if (showBBox)
		renderAABB(fluidBoundingBox, col);
}

void render(const Boundary &boundary, float *boundaryColor)
{
	for (unsigned int i=0; i < boundaries.size(); i++)
	{
		if (renderWalls || (!boundaries[i].isWall))
		{
			meshShader.begin();
			glUniform1f(meshShader.getUniform("shininess"), 5.0f);
			glUniform1f(meshShader.getUniform("specular_factor"), 0.2f);

			GLfloat matrix[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
			glUniformMatrix4fv(meshShader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
			GLfloat pmatrix[16];
			glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
			glUniformMatrix4fv(meshShader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

			glUniform3fv(meshShader.getUniform("surface_color"), 1, &boundaries[i].color[0]);

			MiniGL::drawMesh(boundaries[i].mesh, &boundaries[i].color[0]);

			meshShader.end();
		}
	}
}

void pointShaderBegin(Shader *shader, const float *col, const Real minVal, const Real maxVal, const bool useTexture, float const* color_map)
{
	shader->begin();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glUniform1f(shader->getUniform("viewport_width"), (float)viewport[2]);
	glUniform1f(shader->getUniform("radius"), (float)particleRadius);
	glUniform1f(shader->getUniform("min_scalar"), (GLfloat)minVal);
	glUniform1f(shader->getUniform("max_scalar"), (GLfloat)maxVal);
	glUniform3fv(shader->getUniform("color"), 1, col);

	if (useTexture)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_1D, textureMap);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 256u, 0, GL_RGB, GL_FLOAT, color_map);

		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glGenerateMipmap(GL_TEXTURE_1D);
	}


	GLfloat matrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
	glUniformMatrix4fv(shader->getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
	GLfloat pmatrix[16];
	glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
	glUniformMatrix4fv(shader->getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

	glEnable(GL_DEPTH_TEST);
	// Point sprites do not have to be explicitly enabled since OpenGL 3.2 where
	// they are enabled by default. Moreover GL_POINT_SPRITE is deprecate and only
	// supported before OpenGL 3.2 or with compatibility profile enabled.
	glEnable(GL_POINT_SPRITE);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glPointParameterf(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);
}

void pointShaderEnd(Shader *shader, const bool useTexture)
{
	glBindTexture(GL_TEXTURE_1D, 0);
	shader->end();
}

void updateBoundingBox()
{
	fluidBoundingBox.setEmpty();
	for (size_t i = 0; i < fluids.size(); i++)
	{
		//	const Real r2 = particleRadius*0.5;
		Partio::ParticleAttribute posAttr;
		fluids[i].partioData->attributeInfo(fluids[i].posIndex, posAttr);
		const float* partioX = fluids[i].partioData->data<float>(posAttr, 0);
		for (int j = 0; j < fluids[i].partioData->numParticles(); j++)
		{
			const Eigen::Map<const Eigen::Vector3f> vec(&partioX[3 * j]);
			fluidBoundingBox.extend(vec);
		}
	}
}

void selection(const Vector2i &start, const Vector2i &end, void *clientData)
{
	for (size_t i = 0; i < fluids.size(); i++)
	{
		size_t nParticles;
		if (usePlane)
			nParticles = fluids[i].visibleParticles.size();
		else 
			nParticles = fluids[i].partioData->numParticles();

		if (nParticles != 0)
		{
			Partio::ParticleAttribute posAttr;
			fluids[i].partioData->attributeInfo(fluids[i].posIndex, posAttr);
			const float* partioX = fluids[i].partioData->data<float>(posAttr, 0);
			std::vector<Vector3r> x;
			x.resize(nParticles);
			for (int j = 0; j < nParticles; j++)
			{
				if (usePlane)
					x[j] = Eigen::Map<const Eigen::Vector3f>(&partioX[3 * fluids[i].visibleParticles[j]]).cast<Real>();
				else
					x[j] = Eigen::Map<const Eigen::Vector3f>(&partioX[3 * j]).cast<Real>();
			}

			std::vector<unsigned int> hits;
			fluids[i].selectedParticles.clear();
			Selection::selectRect(start, end, x.begin(),
				x.end(),
				fluids[i].selectedParticles);

			for (int j = 0; j < fluids[i].selectedParticles.size(); j++)
			{
				if (usePlane)
					fluids[i].selectedParticles[j] = fluids[i].visibleParticles[fluids[i].selectedParticles[j]];
				else
					fluids[i].selectedParticles[j] = fluids[i].selectedParticles[j];
			}
		}
	}
}

void particleInfo()
{
	for (size_t i = 0; i < fluids.size(); i++)
	{
		for (unsigned int k = 0; k < fluids[i].selectedParticles.size(); k++)
		{
			unsigned int index = fluids[i].selectedParticles[k];
			std::cout << "Index: " << index << "\n";

			for (int j = 0; j < fluids[i].partioData->numAttributes(); j++)
			{
				Partio::ParticleAttribute attr;
				fluids[i].partioData->attributeInfo(j, attr);

				if (attr.type == Partio::FLOAT)
					std::cout << attr.name << ": " << *fluids[i].partioData->data<float>(attr, index) << "\n";
				else if (attr.type == Partio::INT)
					std::cout << attr.name << ": " << *fluids[i].partioData->data<int>(attr, index) << "\n";
				else if (attr.type == Partio::VECTOR)
				{
					const float *v = fluids[i].partioData->data<float>(attr, index);
					Vector3r vec(v[0], v[1], v[2]);
					std::cout << attr.name << ": " << vec.transpose() << "\n";
				}
			}
			std::cout << "\n";
		}
	}
}

bool readRigidBodyData(std::string fileName, const bool first)
{
	BinaryFileReader binReader;
	if (!binReader.openFile(fileName))
		return false;

	if (first)
	{
		firstRBIndex = frameIndex;
		numberOfRigidBodies = 0;
		binReader.read(numberOfRigidBodies);
		boundaries.resize(numberOfRigidBodies);
		for (unsigned int i = 0; i < numberOfRigidBodies; i++)
		{
			std::string meshFileName = "";
			binReader.read(meshFileName);

			Eigen::Vector3f scale;
			binReader.readMatrix(scale);

			meshFileName = FileSystem::normalizePath(FileSystem::getFilePath(fileName) + "/" + meshFileName);

			loadObj(meshFileName, boundaries[i].mesh, scale.cast<Real>());
			boundaries[i].x0.resize(boundaries[i].mesh.numVertices());
			for (unsigned int j = 0; j < boundaries[i].x0.size(); j++)
				boundaries[i].x0[j] = boundaries[i].mesh.getVertices()[j];

			binReader.read(boundaries[i].isWall);
			binReader.readMatrix(boundaries[i].color);
		}
	}

	for (unsigned int i = 0; i < numberOfRigidBodies; i++)
	{
		binReader.readMatrix(boundaries[i].t);
		binReader.readMatrix(boundaries[i].R);
 		for (unsigned int j = 0; j < boundaries[i].x0.size(); j++)
			boundaries[i].mesh.getVertices()[j] = boundaries[i].R.transpose().cast<Real>() * boundaries[i].x0[j] + boundaries[i].t.cast<Real>();

		boundaries[i].mesh.updateNormals();
		boundaries[i].mesh.updateVertexNormals();
	}

	binReader.closeFile();
	return true;
}


bool readPartioFile(const std::string &fileName, Partio::ParticlesDataMutable* &partioData, unsigned int &posIndex, const bool printInfo)
{
	if (!FileSystem::fileExists(fileName))
		return false;

	partioData = Partio::read(fileName.c_str());

	if (!partioData)
		return false;

	LOG_INFO << "Frame: " << frameIndex << ", number of particles: " << partioData->numParticles();

	posIndex = 0xffffffff;

	for (int i = 0; i < partioData->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		partioData->attributeInfo(i, attr);
		if (attr.name == "position")
			posIndex = i;
		if (printInfo)
			LOG_INFO << "Found attribute: " << attr.name;
	}
	if (posIndex == 0xffffffff)
		return false;

	return true;
}


std::string zeroPadding(const unsigned int number, const unsigned int length)
{
	std::ostringstream out;
	out << std::internal << std::setfill('0') << std::setw(length) << number;
	return out.str();
}

unsigned int getFrameIndexFromFile(const std::string &inputFileName)
{
	std::string fileName = inputFileName;

	std::smatch match;
	std::string tmp = fileName;
	std::string frameNumber = "";
	while (regex_search(tmp, match, std::regex("([0-9]+)")))
	{
		frameNumber = match.str();
		tmp = match.suffix();
	}
	if (frameNumber != "")
		return (unsigned int) std::stoi(frameNumber);
	else
		return 0xffffffff;
}

std::string convertFileName(const std::string &inputFileName, const std::string currentFrame)
{
	std::string fileName = inputFileName;

	std::smatch match;
	std::string tmp = fileName;
	std::string frameNumber = "";
	while (regex_search(tmp, match, std::regex("([0-9]+)")))
	{
		frameNumber = match.str();
		tmp = match.suffix();
	}
	std::string res = "";
	if (frameNumber != "")
		res = frameNumber;
	else
		return "";

	std::string::size_type pos1 = fileName.find_last_of(res) - res.size() + 1;

	//std::string numberStr = zeroPadding(currentFrame, (unsigned int)res.size());
	std::string numberStr = currentFrame;
	fileName.replace(pos1, res.size(), numberStr);
	return fileName;
}

bool nextFrame()
{
	frameIndex++;
	if (!updateData())
	{
		LOG_INFO << "File with index " << frameIndex << " does not exist.";
		frameIndex--;
		return false;
	}
	TwRefreshBar(MiniGL::getTweakBar());
	return true;
}

bool prevFrame()
{
	if (frameIndex > 0)
	{
		frameIndex--;
		if (!updateData())
		{
			LOG_INFO << "File with index " << frameIndex << " does not exist.";
			frameIndex++;
			return false;
		}
	}
	TwRefreshBar(MiniGL::getTweakBar());
	return true;
}


bool setFrame(const unsigned int index)
{
	unsigned int lastIndex = frameIndex;
	frameIndex = index;
	if (!updateData())
	{
		LOG_INFO << "File with index " << frameIndex << " does not exist.";
		frameIndex = lastIndex;
		return false;
	}
	TwRefreshBar(MiniGL::getTweakBar());
	return true;
}

bool updateData()
{
	bool chk = false;
	for (size_t i = 0; i < fluids.size(); i++)
	{
		fluids[i].currentFile = convertFileName(fluids[i].inputFile, to_string(frameIndex));
		LOG_INFO << fluids[i].currentFile;
		if (readPartioFile(fluids[i].currentFile, fluids[i].partioData, fluids[i].posIndex))
			chk = true;

		if (showBBox)
			updateBoundingBox();

		if (usePlane)
		{
			Partio::ParticleAttribute posAttr;
			fluids[i].partioData->attributeInfo(fluids[i].posIndex, posAttr);
			const float* partioX = fluids[i].partioData->data<float>(posAttr, 0);
			fluids[i].visibleParticles.clear();
			fluids[i].visibleParticles.reserve(fluids[i].partioData->numParticles());
			Vector3f normal = planeNormal;
			normal.normalize();
			for (unsigned int j = 0; j < fluids[i].partioData->numParticles(); j++)
			{
				const Eigen::Map<const Eigen::Vector3f> vec(&partioX[3 * j]);
				if ((vec.dot(normal) - planePoint.dot(normal)) > 0)
					fluids[i].visibleParticles.push_back(j);
			}
		}
	}
	if (chk)
	{
		if (useRBData)
		{
			string currentFile = convertFileName(rbDataFile, to_string(frameIndex));
			LOG_INFO << currentFile;
			readRigidBodyData(currentFile, frameIndex == firstRBIndex);
		}
	}
	return chk;
}

bool imagesExist(const unsigned int frameIndex)
{
	bool res = true;
	for (size_t i = 0; i < fluids.size(); i++)
	{
		const std::string currentFile = convertFileName(fluids[i].inputFile, to_string(frameIndex));
		const std::string imageFileName = FileSystem::normalizePath(outPath + "/" + FileSystem::getFileName(currentFile) + ".jpg");
		if (!FileSystem::fileExists(imageFileName))
			res = false;
	}
	return res;
}

void getImage()
{
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadBuffer(GL_FRONT_LEFT);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, image.data());
	flipImage(width, height, image.data());
}

void addImageDataToVideo()
{
	fwrite(image.data(), width*height * 3, 1, ffmpegPipe);
}

void saveImageData()
{
	FileSystem::makeDirs(outPath);
	std::string fileName = FileSystem::normalizePath(outPath + "/" + FileSystem::getFileName(fluids[0].currentFile) + ".jpg");
	LOG_INFO << "Write file: " << fileName;
	saveImage(fileName);
}

void saveFrame()
{
	MiniGL::viewport();
	render();
	glutSwapBuffers();
	getImage();
	saveImageData();
}

void outputJpg(unsigned char oneByte) 
{ 
	fputc(oneByte, jpegFile);
}

void flipImage(int width, int height, unsigned char *image)
{
	unsigned char tmp[3];

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height / 2; y++)
		{
			int top = (x + y * width) * 3;
			int bottom = (x + (height - y - 1) * width) * 3;

			memcpy(tmp, image + top, sizeof(tmp));
			memcpy(image + top, image + bottom, sizeof(tmp));
			memcpy(image + bottom, tmp, sizeof(tmp));
		}
	}
}

void saveImage(const std::string &fileName)
{
 	jpegFile = fopen(fileName.c_str(), "wb");
 	TooJpeg::writeJpeg(outputJpg, image.data(), width, height);
 	fclose(jpegFile);
}

void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale)
{
	std::vector<OBJLoader::Vec3f> x;
	std::vector<OBJLoader::Vec3f> normals;
	std::vector<MeshFaceIndices> faces;
	OBJLoader::Vec3f s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	OBJLoader::loadObj(filename, &x, &faces, &normals, nullptr, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	mesh.initMesh(nPoints, nFaces);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		mesh.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		// Reduce the indices by one
		int posIndices[3];
		for (int j = 0; j < 3; j++)
		{
			posIndices[j] = faces[i].posIndices[j] - 1;
		}

		mesh.addFace(&posIndices[0]);
	}

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;
}


void TW_CALL setFrameIndex(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	setFrame(val);
}

void TW_CALL getFrameIndex(void *value, void *clientData)
{
	*(unsigned int *)(value) = frameIndex;
}

void TW_CALL setShowBBox(const void *value, void *clientData)
{
	const bool val = *(const bool *)(value);
	showBBox = val;
	if (showBBox)
		updateBoundingBox();
}

void TW_CALL getShowBBox(void *value, void *clientData)
{
	*(bool *)(value) = showBBox;
}

void TW_CALL setUsePlane(const void *value, void *clientData)
{
	const bool val = *(const bool *)(value);
	usePlane = val;
	if (usePlane)
		updateData();
}

void TW_CALL getUsePlane(void *value, void *clientData)
{
	*(bool *)(value) = usePlane;
}

void TW_CALL setPlaneNormal(const void *value, void *clientData)
{
	float * const val = (float* const)(value);
	planeNormal = Vector3f(val[0], val[1], val[2]);
	if (usePlane)
		updateData();
}

void TW_CALL getPlaneNormal(void *value, void *clientData)
{
	((float*)value)[0] = planeNormal[0];
	((float*)value)[1] = planeNormal[1];
	((float*)value)[2] = planeNormal[2];
}

void TW_CALL setPlanePoint(const void *value, void *clientData)
{
	float * const val = (float* const)(value);
	planePoint = Vector3f(val[0], val[1], val[2]);
	if (usePlane)
		updateData();
}

void TW_CALL getPlanePoint(void *value, void *clientData)
{
	((float*)value)[0] = planePoint[0];
	((float*)value)[1] = planePoint[1];
	((float*)value)[2] = planePoint[2];
}
