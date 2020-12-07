#include "PartioViewer.h"

#include "Utilities/OBJLoader.h"
#include "Utilities/FileSystem.h"
#include "Utilities/BinaryFileReaderWriter.h"
#include "Utilities/Version.h"
#include "extern/cxxopts/cxxopts.hpp"
#include "Utilities/Timing.h"
#include "extern/toojpeg/toojpeg.h"
#include "GUI/OpenGL/PartioViewer_OpenGL.h"
#ifdef USE_IMGUI
#include "GUI/imgui/PartioViewer_GUI_imgui.h"
#else
#include "GUI/TweakBar/PartioViewer_GUI_TweakBar.h"
#endif

using namespace std;
using namespace SPH;
using namespace Utilities;

FILE *PartioViewer::m_jpegFile = nullptr;


PartioViewer::PartioViewer()
{
	m_frameIndex = 1;
	m_startFrame = -1;
	m_endFrame = -1;
	m_firstRBIndex = 1;
	m_numberOfRigidBodies = 0;
	m_particleRadius = 0.025;
	m_planeNormal = Vector3f(0, 0, -1);
	m_planePoint.setZero();
	m_renderSequence = false;
	m_renderVideo = false;
	m_overWrite = true;
	m_ffmpegPath = "ffmpeg";
	m_rbDataFile = "";
	m_ffmpegPipe = nullptr;
	m_width = 1280;
	m_height = 960;
	m_fps = 25;
	m_doPause = true;
	m_usePlane = false;
	m_useRBData = false;
	m_jpegFile = nullptr;
}


PartioViewer::~PartioViewer()
{
}

int PartioViewer::run(int argc, char **argv)
{
	m_argc = argc;
	m_argv = argv;

	Utilities::logger.addSink(unique_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));

	LOG_INFO << "Git refspec: " << GIT_REFSPEC;
	LOG_INFO << "Git SHA1: " << GIT_SHA1;
	LOG_INFO << "Git status: " << GIT_LOCAL_STATUS;

	m_exePath = FileSystem::getProgramPath();
#ifdef USE_IMGUI
	m_gui = new PartioViewer_GUI_imgui(this);
#else
	m_gui = new PartioViewer_GUI_TweakBar(this);
#endif

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
			("colorField", "Name of field that is used as default for the color.", cxxopts::value<std::string>()->default_value("velocity"))
			("colorMapType", "Default color map (0=None, 1=Jet, 2=Plasma, 3=CoolWarm, 4=BlueWhiteRed, 5=Seismic)", cxxopts::value<unsigned int>()->default_value("1"))
			("renderMinValue", "Default min value of field.", cxxopts::value<float>()->default_value("0.0"))
			("renderMaxValue", "Default max value of field.", cxxopts::value<float>()->default_value("10.0"))
			;

		m_gui->addOptions(options);

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
								m_rbDataFile = FileSystem::normalizePath(rbDir + "/" + framePairs[0].second);
								m_useRBData = true;
							}
						}
					}
				}
			}


			getFluids().resize(inputFile.size());
			for (int i = 0; i < inputFile.size(); i++)
			{
				getFluids()[i].inputFile = inputFile[i];
				if (FileSystem::isRelativePath(inputFile[i]))
					getFluids()[i].inputFile = FileSystem::normalizePath(m_exePath + "/" + inputFile[i]);
				getFluids()[i].currentFile = getFluids()[i].inputFile;
			}
		}
		else
		{
			LOG_INFO << options.help({ "", "Group" });
			exit(0);
		}

		m_outPath = m_exePath;
		if (result.count("outdir"))
			m_outPath = result["outdir"].as<std::string>();

		if (result.count("rbData"))
		{
			m_rbDataFile = result["rbData"].as<std::string>();
			m_useRBData = true;
		}

		if (result.count("ffmpegPath"))
			m_ffmpegPath = result["ffmpegPath"].as<std::string>();

		m_particleRadius = 0.025;
		if (result.count("radius"))
			m_particleRadius = result["radius"].as<Real>();
		LOG_INFO << "Radius: " << m_particleRadius;

		if (result.count("startFrame"))
			m_startFrame = result["startFrame"].as<int>();
		LOG_INFO << "Start frame: " << m_startFrame;

		if (result.count("endFrame"))
			m_endFrame = result["endFrame"].as<int>();
		LOG_INFO << "End frame: " << m_endFrame;

		if (result.count("renderSequence"))
		{
			m_renderSequence = true;
			if (result.count("noOverwrite"))
				m_overWrite = false;
		}

		if (result.count("renderVideo"))
			m_renderVideo = true;

		if (result.count("width"))
			m_width = result["width"].as<int>();

		if (result.count("height"))
			m_height = result["height"].as<int>();

		if (result.count("fps"))
			m_fps = result["fps"].as<int>();

		m_defaultColorFieldName = "velocity";
		if (result.count("colorField"))
			m_defaultColorFieldName = result["colorField"].as<std::string>();

		m_defaultColorMapType = 1;
		if (result.count("colorMapType"))
			m_defaultColorMapType = result["colorMapType"].as<unsigned int>();

		m_defaultRenderMinValue = 0.0;
		if (result.count("renderMinValue"))
			m_defaultRenderMinValue = result["renderMinValue"].as<float>();

		m_defaultRenderMaxValue = 10.0;
		if (result.count("renderMaxValue"))
			m_defaultRenderMaxValue = result["renderMaxValue"].as<float>();

		m_gui->parseOptions(result);

		m_image.resize(3 * m_width * m_height);
	}
	catch (const cxxopts::OptionException& e)
	{
		LOG_INFO << "error parsing options: " << e.what();
		exit(1);
	}

	m_frameIndex = getFrameIndexFromFile(getFluids()[0].inputFile);
	if (m_frameIndex == 0xffffffff)
		m_frameIndex = 0;
	m_firstRBIndex = m_frameIndex;

	// read particle data
	if (!updateData())
	{
		LOG_ERR << "Unable to read partio file or did not find position data.";
		exit(1);
	}

	// set default color field
	for (auto i = 0; i < m_fluids.size(); i++)
	{
		if (m_fluids[i].partioData)
		{
			bool found = false;
			int idx = 0;
			int velIdx = 0;
			for (int j = 0; j < m_fluids[i].partioData->numAttributes(); j++)
			{
				Partio::ParticleAttribute attr;
				m_fluids[i].partioData->attributeInfo(j, attr);
				if ((attr.type == Partio::FLOAT) || (attr.type == Partio::VECTOR) || (attr.type == Partio::INT))
				{
					// select the default color field name
					if (attr.name == m_defaultColorFieldName)
					{
						m_fluids[i].m_colorField = idx;
						found = true;
						break;
					}
					if (attr.name == "velocity")
						velIdx = idx;

					idx++;
				}
			}
			// if default color field name was not found, use velocity as default
			if (!found)
				m_fluids[i].m_colorField = velIdx;
			m_fluids[i].m_colorMapType = m_defaultColorMapType;
			m_fluids[i].m_renderMinValue = m_defaultRenderMinValue;
			m_fluids[i].m_renderMaxValue = m_defaultRenderMaxValue;
		}
	}


	m_gui->init();

	initGUI();
	updateBoundingBox();

	if ((m_renderSequence || m_renderVideo) && (m_startFrame >= 0))
		setFrame(m_startFrame);

	m_gui->run();

	m_gui->cleanup();

	for (size_t i = 0; i < getFluids().size(); i++)
		getFluids()[i].partioData->release();
	Timing::printAverageTimes();
	Timing::printTimeSums();

	delete m_gui;

	return 0;
}

void PartioViewer::timeStep()
{
	if (!m_renderSequence && !m_renderVideo)
	{
		if (!m_doPause)
		{
			if (!nextFrame())
				m_doPause = true;
		}
	}
	else if (m_renderVideo)
	{
		generateVideo();
		m_gui->stop();
	}
	else if (m_renderSequence)
	{
		generateSequence();
		m_gui->stop();
	}
}

void PartioViewer::loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale)
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


bool PartioViewer::readPartioFile(const std::string &fileName, Partio::ParticlesDataMutable* &partioData, unsigned int &posIndex, const bool printInfo)
{
	if (!FileSystem::fileExists(fileName))
		return false;

	if (partioData)
		partioData->release();
	partioData = Partio::read(fileName.c_str());

	if (!partioData)
		return false;

	LOG_INFO << "Frame: " << m_frameIndex << ", number of particles: " << partioData->numParticles();

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

bool PartioViewer::readRigidBodyData(std::string fileName, const bool first)
{
	BinaryFileReader binReader;
	if (!binReader.openFile(fileName))
		return false;

	if (first)
	{
		m_firstRBIndex = m_frameIndex;
		m_numberOfRigidBodies = 0;
		binReader.read(m_numberOfRigidBodies);
		m_boundaries.resize(m_numberOfRigidBodies);
		for (unsigned int i = 0; i < m_numberOfRigidBodies; i++)
		{
			std::string meshFileName = "";
			binReader.read(meshFileName);

			Eigen::Vector3f scale;
			binReader.readMatrix(scale);

			meshFileName = FileSystem::normalizePath(FileSystem::getFilePath(fileName) + "/" + meshFileName);

			loadObj(meshFileName, m_boundaries[i].mesh, scale.cast<Real>());
			m_boundaries[i].x0.resize(m_boundaries[i].mesh.numVertices());
			for (unsigned int j = 0; j < m_boundaries[i].x0.size(); j++)
				m_boundaries[i].x0[j] = m_boundaries[i].mesh.getVertices()[j];

			binReader.read(m_boundaries[i].isWall);
			binReader.readMatrix(m_boundaries[i].color);
		}
	}

	for (unsigned int i = 0; i < m_numberOfRigidBodies; i++)
	{
		binReader.readMatrix(m_boundaries[i].t);
		binReader.readMatrix(m_boundaries[i].R);
		for (unsigned int j = 0; j < m_boundaries[i].x0.size(); j++)
			m_boundaries[i].mesh.getVertices()[j] = m_boundaries[i].R.cast<Real>() * m_boundaries[i].x0[j] + m_boundaries[i].t.cast<Real>();

		m_boundaries[i].mesh.updateNormals();
		m_boundaries[i].mesh.updateVertexNormals();
	}

	binReader.closeFile();
	return true;
}


unsigned int PartioViewer::getFrameIndexFromFile(const std::string &inputFileName)
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
		return (unsigned int)std::stoi(frameNumber);
	else
		return 0xffffffff;
}

std::string PartioViewer::convertFileName(const std::string &inputFileName, const std::string currentFrame)
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


bool PartioViewer::nextFrame()
{
	m_frameIndex++;
	if (!updateData())
	{
		LOG_INFO << "File with index " << m_frameIndex << " does not exist.";
		m_frameIndex--;
		return false;
	}
	m_gui->update();
	return true;
}

bool PartioViewer::prevFrame()
{
	if (m_frameIndex > 0)
	{
		m_frameIndex--;
		if (!updateData())
		{
			LOG_INFO << "File with index " << m_frameIndex << " does not exist.";
			m_frameIndex++;
			return false;
		}
	}
	m_gui->update();
	return true;
}


bool PartioViewer::setFrame(const unsigned int index)
{
	unsigned int lastIndex = m_frameIndex;
	m_frameIndex = index;
	if (!updateData())
	{
		LOG_INFO << "File with index " << m_frameIndex << " does not exist.";
		m_frameIndex = lastIndex;
		return false;
	}
	m_gui->update();
	return true;
}

std::string PartioViewer::zeroPadding(const unsigned int number, const unsigned int length)
{
	std::ostringstream out;
	out << std::internal << std::setfill('0') << std::setw(length) << number;
	return out.str();
}


bool PartioViewer::updateData()
{
	bool chk = false;
	for (size_t i = 0; i < getFluids().size(); i++)
	{
		getFluids()[i].currentFile = convertFileName(getFluids()[i].inputFile, to_string(m_frameIndex));
		LOG_INFO << getFluids()[i].currentFile;
		if (readPartioFile(getFluids()[i].currentFile, getFluids()[i].partioData, getFluids()[i].posIndex))
			chk = true;

		updateBoundingBox();

		if (m_usePlane)
		{
			Partio::ParticleAttribute posAttr;
			getFluids()[i].partioData->attributeInfo(getFluids()[i].posIndex, posAttr);
			const float* partioX = getFluids()[i].partioData->data<float>(posAttr, 0);
			getFluids()[i].visibleParticles.clear();
			getFluids()[i].visibleParticles.reserve(getFluids()[i].partioData->numParticles());
			Vector3f normal = m_planeNormal;
			normal.normalize();
			for (auto j = 0; j < getFluids()[i].partioData->numParticles(); j++)
			{
				const Eigen::Map<const Eigen::Vector3f> vec(&partioX[3 * j]);
				if ((vec.dot(normal) - m_planePoint.dot(normal)) > 0)
					getFluids()[i].visibleParticles.push_back(j);
			}
		}
	}
	if (chk)
	{
		if (m_useRBData)
		{
			string currentFile = convertFileName(m_rbDataFile, to_string(m_frameIndex));
			LOG_INFO << currentFile;
			readRigidBodyData(currentFile, m_frameIndex == m_firstRBIndex);
		}
		PartioViewer_OpenGL::updateScalarField();
	}
	return chk;
}

bool PartioViewer::imagesExist(const unsigned int frameIndex)
{
	bool res = true;
	for (size_t i = 0; i < getFluids().size(); i++)
	{
		const std::string currentFile = convertFileName(getFluids()[i].inputFile, to_string(m_frameIndex));
		const std::string imageFileName = FileSystem::normalizePath(m_outPath + "/" + FileSystem::getFileName(currentFile) + ".jpg");
		if (!FileSystem::fileExists(imageFileName))
			res = false;
	}
	return res;
}

void PartioViewer::addImageDataToVideo()
{
	fwrite(m_image.data(), m_width*m_height * 3, 1, m_ffmpegPipe);
}

void PartioViewer::saveImageData()
{
	FileSystem::makeDirs(m_outPath);
	std::string fileName = FileSystem::normalizePath(m_outPath + "/" + FileSystem::getFileName(getFluids()[0].currentFile) + ".jpg");
	LOG_INFO << "Write file: " << fileName;
	saveImage(fileName, m_image.data(), m_width, m_height);
}

void PartioViewer::saveFrame()
{
	m_gui->render();
	PartioViewer_OpenGL::getImage(m_width, m_height, m_image.data());
	saveImageData();
}

void PartioViewer::outputJpg(unsigned char oneByte)
{
	fputc(oneByte, m_jpegFile);
}


void PartioViewer::saveImage(const std::string &fileName, const void *pixels, unsigned int width, unsigned int height)
{
	m_jpegFile = fopen(fileName.c_str(), "wb");
	TooJpeg::writeJpeg(outputJpg, pixels, width, height);
	fclose(m_jpegFile);
}


void PartioViewer::generateSequence()
{
	m_width = m_gui->getWidth();
	m_height = m_gui->getHeight();

	// height must be divisible by 2
	if (m_height % 2 == 1)
		m_height--;
	// width must be divisible by 2
	if (m_width % 2 == 1)
		m_width--;

	m_image.resize(3 * m_width*m_height);

	if ((m_startFrame >= 0) && (!imagesExist(m_startFrame)))
		setFrame(m_startFrame);

	bool chk = true;
	while (chk && ((m_frameIndex <= m_endFrame) || (m_endFrame < 0)))
	{
		if (m_overWrite || !imagesExist(m_frameIndex))
		{
			if (!updateData())
			{
				m_frameIndex--;
				chk = false;
				break;
			}

			m_gui->render();

			PartioViewer_OpenGL::getImage(m_width, m_height, m_image.data());
			saveImageData();
		}
		m_frameIndex++;
	}
}


void PartioViewer::generateVideo()
{
	m_width = m_gui->getWidth();
	m_height = m_gui->getHeight();

	// height must be divisible by 2
	if (m_height % 2 == 1)
		m_height--;
	// width must be divisible by 2
	if (m_width % 2 == 1)
		m_width--;

	// open pipe to ffmpeg
	stringstream sstm;

	FileSystem::makeDirs(m_outPath);
	std::string videoFile = FileSystem::normalizePath(m_outPath + "/output.mp4");
	sstm << m_ffmpegPath << " -y -hide_banner -nostats -loglevel panic -r " << m_fps << " -f rawvideo -vcodec rawvideo -pix_fmt rgb24 -s "
		<< m_width << "x" << m_height << " -i - -threads 0 -preset fast -pix_fmt yuv420p -crf 21 " << videoFile << " > ffmpeg.log";

	#ifdef WIN32
	if (!(m_ffmpegPipe = _popen(sstm.str().c_str(), "wb")))
		#else
	if (!(m_ffmpegPipe = popen(sstm.str().c_str(), "w")))
		#endif
	{
		LOG_ERR << "Cannot open pipe to ffmpeg.";
		return;
	}

	m_image.resize(3 * m_width*m_height);

	if (m_startFrame >= 0)
		setFrame(m_startFrame);

	bool chk = true;
	while (chk && ((m_frameIndex <= m_endFrame) || (m_endFrame < 0)))
	{
		m_gui->render();

		PartioViewer_OpenGL::getImage(m_width, m_height, m_image.data());
		addImageDataToVideo();
		chk = nextFrame();
	}
	fflush(m_ffmpegPipe);
	fclose(m_ffmpegPipe);
	m_ffmpegPipe = nullptr;
	LOG_INFO << "Generated video: " << videoFile;
}


void PartioViewer::particleInfo()
{
	for (size_t i = 0; i < getFluids().size(); i++)
	{
		for (unsigned int k = 0; k < getFluids()[i].selectedParticles.size(); k++)
		{
			unsigned int index = getFluids()[i].selectedParticles[k];
			std::cout << "Index: " << index << "\n";

			for (int j = 0; j < getFluids()[i].partioData->numAttributes(); j++)
			{
				Partio::ParticleAttribute attr;
				getFluids()[i].partioData->attributeInfo(j, attr);

				if (attr.type == Partio::FLOAT)
					std::cout << attr.name << ": " << *getFluids()[i].partioData->data<float>(attr, index) << "\n";
				else if (attr.type == Partio::INT)
					std::cout << attr.name << ": " << *getFluids()[i].partioData->data<int>(attr, index) << "\n";
				else if (attr.type == Partio::VECTOR)
				{
					const float *v = getFluids()[i].partioData->data<float>(attr, index);
					Vector3r vec(v[0], v[1], v[2]);
					std::cout << attr.name << ": " << vec.transpose() << "\n";
				}
			}
			std::cout << "\n";
		}
	}
}

void PartioViewer::updateBoundingBox()
{
	m_fluidBoundingBox.setEmpty();
	for (size_t i = 0; i < m_fluids.size(); i++)
	{
		if (m_fluids[i].partioData)
		{
			//	const Real r2 = particleRadius*0.5;
			Partio::ParticleAttribute posAttr;
			m_fluids[i].partioData->attributeInfo(m_fluids[i].posIndex, posAttr);
			const float* partioX = m_fluids[i].partioData->data<float>(posAttr, 0);
			for (int j = 0; j < m_fluids[i].partioData->numParticles(); j++)
			{
				const Eigen::Map<const Eigen::Vector3f> vec(&partioX[3 * j]);
				if (!std::isfinite(vec[0]) || !std::isfinite(vec[1]) || !std::isfinite(vec[2]))
					LOG_WARN << "Warning: Fluid " << i << ", Particle " << j << ": " << vec.transpose();
				else
					m_fluidBoundingBox.extend(vec);
			}
		}
	}
	m_fluidBoundingBox.extend(m_fluidBoundingBox.max() + m_particleRadius * Vector3f::Ones());
	m_fluidBoundingBox.extend(m_fluidBoundingBox.min() - m_particleRadius * Vector3f::Ones());
}

void PartioViewer::initGUI()
{
	m_gui->initParameterGUI();
}

void PartioViewer::setPlanePoint(const Eigen::Vector3f &val)
{
	m_planePoint = val; 
	if (m_usePlane)
		updateData();
}

void PartioViewer::setPlaneNormal(const Eigen::Vector3f &val)
{
	m_planeNormal = val;
	if (m_usePlane)
		updateData();
}

void PartioViewer::setUsePlane(bool val)
{
	m_usePlane = val;
	if (m_usePlane)
		updateData();
}

void SPH::PartioViewer::reset()
{
	setFrame(m_firstRBIndex);
}
