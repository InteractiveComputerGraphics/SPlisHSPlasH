#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/FileSystem.h"
#include "Utilities/Version.h"
#include "extern/partio/src/lib/Partio.h"
#include "extern/cxxopts/cxxopts.hpp"

#include <regex>

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

INIT_TIMING
INIT_LOGGING

using namespace Utilities;

void saveParticleCloudVTK(const std::string & path, const Partio::ParticlesDataMutable * partioData);

std::string inputFile = "";
int startFrame = -1;
int endFrame = -1;
std::string exePath, dataPath, outDir;

template<typename T>
inline void swapByteOrder(T*v)
{
	constexpr size_t n = sizeof(T);
	uint8_t out[n];
	for (unsigned int c = 0; c < n; c++)
		out[c] = reinterpret_cast<uint8_t*>(v)[n - c - 1];
	std::memcpy(v, out, n);
}

// main 
int main(int argc, char **argv)
{
	REPORT_MEMORY_LEAKS;

	Utilities::logger.addSink(std::unique_ptr<Utilities::ConsoleSink>(new Utilities::ConsoleSink(Utilities::LogLevel::INFO)));

	LOG_INFO << "Git refspec: " << GIT_REFSPEC;
	LOG_INFO << "Git SHA1: " << GIT_SHA1;
	LOG_INFO << "Git status: " << GIT_LOCAL_STATUS;

	exePath = FileSystem::getProgramPath();
	dataPath = FileSystem::normalizePath(exePath + "/" + std::string(SPH_DATA_PATH));

	try
	{
		cxxopts::Options options(argv[0], "partio2vtk - Converts partio to vtk files");
		options
			.positional_help("[single-file or sequence, e.g. particles_#.bgeo]")
			.show_positional_help();

		options.add_options()
			("h,help", "Print help")
			("s,startFrame", "Start frame (only used if value is >= 0)", cxxopts::value<int>()->default_value("-1"))
			("e,endFrame", "End frame (only used if value is >= 0)", cxxopts::value<int>()->default_value("-1"))
			("o,outputDir", "Output directory", cxxopts::value<std::string>())
			;

		options.add_options("invisible")
			("input-file", "Input file", cxxopts::value<std::string>());

		options.parse_positional("input-file");
		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "" }) << std::endl;
			exit(0);
		}

		if (result.count("input-file"))
		{
			inputFile = result["input-file"].as<std::string>();;
		}
		else
		{
			std::cout << options.help({ "" }) << std::endl;
			exit(0);
		}

		if (result.count("outputDir"))
			outDir = result["outputDir"].as<std::string>();
		LOG_INFO << "Output directory: " << outDir;

		if (result.count("startFrame"))
			startFrame = result["startFrame"].as<int>();
		LOG_INFO << "Start frame: " << startFrame;

		if (result.count("endFrame"))
			endFrame = result["endFrame"].as<int>();
		LOG_INFO << "End frame: " << endFrame;
		
	}
	catch (const cxxopts::OptionException& e)
	{
		LOG_INFO << "error parsing options: " << e.what();
		exit(1);
	}
	
	const std::string fileName = FileSystem::getFileName(inputFile);
	const std::string filePath = FileSystem::getFilePath(inputFile);
	const std::string fileExt = FileSystem::getFileExt(inputFile);

	FileSystem::makeDirs(outDir);

	const auto hashPos = fileName.find_first_of('#');
	// convert single file
	if (std::string::npos == hashPos)
	{
		LOG_INFO << "Converting a single file";
		Partio::ParticlesDataMutable * partioData = Partio::read(inputFile.c_str());
		if (nullptr == partioData)
		{
			LOG_ERR << "Could not read file " << inputFile;
			return -1;
		}
		LOG_INFO << "Successfully read file " << inputFile;
		std::string outputFile = outDir + "/" + fileName + ".vtk";
		LOG_INFO << "Writing file " << outputFile;
		saveParticleCloudVTK(outputFile, partioData);
		LOG_INFO << "Freeing data";
		partioData->release();
		LOG_INFO << "Done";
	}
	// convert files that match the pattern
	else
	{
		LOG_INFO << "Converting a range of files";
		const auto hashesEnd = fileName.find_last_of('#');
		unsigned int start = 1;
		if (startFrame >= 0)
			start = startFrame;
		unsigned int end = std::numeric_limits<unsigned int>::max();
		if (endFrame >= 0)
			end = endFrame;
		for (unsigned int i = start; i <= end; i++)
		{
			std::string fileNameWithNumber;
			std::regex_replace(std::back_inserter(fileNameWithNumber), fileName.begin(), fileName.end(), std::regex("#+"), std::to_string(i));
			std::string inputFileWithNumber = filePath + "/" + fileNameWithNumber + "." + fileExt;

			if (!FileSystem::fileExists(inputFileWithNumber))
			{
				LOG_INFO << "File " << inputFileWithNumber << " does not exist. Assuming the last file has been read.";
				break;
			}
			Partio::ParticlesDataMutable * partioData = Partio::read(inputFileWithNumber.c_str());
			if (nullptr == partioData)
				break;

			LOG_INFO << "Successfully read file " << inputFileWithNumber;
			std::string outputFile = outDir + "/" + fileNameWithNumber + ".vtk";
			LOG_INFO << "Writing file " << outputFile;
			saveParticleCloudVTK(outputFile, partioData);
			LOG_INFO << "Freeing data";
			partioData->release();
		}
	}

	Timing::printAverageTimes();
	Timing::printTimeSums();
	
	return 0;
}

/*
Note: Binary VTK works with big endianness.
*/
void saveParticleCloudVTK(const std::string & path, const Partio::ParticlesDataMutable * partioData)
{
	const unsigned int numParticles = partioData->numParticles();
	if (0 == numParticles)
		return;

	START_TIMING("Writing VTK file");
	// Open the file
	std::ofstream outfile{ path, std::ios::binary };
	if (!outfile.is_open()) {
		LOG_ERR << "Cannot open a file to save a VTK mesh.";
		exit(-1);
	}

	outfile << "# vtk DataFile Version 4.1\n";
	outfile << "\n";
	outfile << "BINARY\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";

	//////////////////////////////////////////////////////////////////////////
	// find indices of position and ID attribute
	unsigned int posIndex = 0xffffffff;
	unsigned int idIndex = 0xffffffff;
	for (int i = 0; i < partioData->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		partioData->attributeInfo(i, attr);
		if (attr.name == "position")
			posIndex = i;
		else if (attr.name == "id")
			idIndex = i;

		LOG_INFO << "Found attribute: " << attr.name;
	}

	//////////////////////////////////////////////////////////////////////////
	// export position attribute as POINTS
	if (0xffffffff != posIndex)
	{
		// copy from partio data
		std::vector<Vector3f> positions;
		positions.reserve(numParticles);
		Partio::ParticleAttribute attr;
		partioData->attributeInfo(posIndex, attr);
		for (unsigned int i = 0u; i < numParticles; i++)
			positions.emplace_back(partioData->data<float>(attr, i));
		// swap endianess
		for (unsigned int i = 0; i < numParticles; i++)
			for (unsigned int c = 0; c < 3; c++)
				swapByteOrder(&positions[i][c]);
		// export to vtk
		outfile << "POINTS " << numParticles << " float\n";
		outfile.write(reinterpret_cast<char*>(positions[0].data()), 3 * numParticles * sizeof(float));
		outfile << "\n";
	}
	else
	{
		LOG_ERR << "No particle positions found!";
		STOP_TIMING_AVG;
		return;
	}

	//////////////////////////////////////////////////////////////////////////
	// export particle IDs as CELLS
	{
		std::vector<Eigen::Vector2i> cells;
		cells.reserve(numParticles);
		int nodes_per_cell_swapped = 1;
		swapByteOrder(&nodes_per_cell_swapped);
		if (0xffffffff != idIndex)
		{
			// load IDs from partio
			Partio::ParticleAttribute attr;
			partioData->attributeInfo(idIndex, attr);
			for (unsigned int i = 0u; i < numParticles; i++)
			{
				int idSwapped = *partioData->data<int>(attr, i);
				swapByteOrder(&idSwapped);
				cells.emplace_back(nodes_per_cell_swapped, idSwapped);
			}
		}
		else
		{
			// generate IDs
			for (unsigned int i = 0u; i < numParticles; i++)
			{
				int idSwapped = i;
				swapByteOrder(&idSwapped);
				cells.emplace_back(nodes_per_cell_swapped, idSwapped);
			}
		}

		// particles are cells with one element and the index of the particle
		outfile << "CELLS " << numParticles << " " << 2 * numParticles << "\n";
		outfile.write(reinterpret_cast<char*>(cells[0].data()), 2 * numParticles * sizeof(int));
		outfile << "\n";
	}
	//////////////////////////////////////////////////////////////////////////
	// export cell types
	{
		// the type of a particle cell is always 1
		std::vector<int> cellTypes;
		unsigned int cellTypeSwapped = 1;
		swapByteOrder(&cellTypeSwapped);
		cellTypes.resize(numParticles, cellTypeSwapped);
		outfile << "CELL_TYPES " << numParticles << "\n";
		outfile.write(reinterpret_cast<char*>(cellTypes.data()), numParticles * sizeof(int));
		outfile << "\n";
	}
	//////////////////////////////////////////////////////////////////////////
	// write additional attributes as per-particle data
	outfile << "POINT_DATA " << numParticles << "\n";
	// per point fields (all attributes except for positions and IDs)
	const unsigned int numFields = partioData->numAttributes() - static_cast<int>(0xffffffff != posIndex) - static_cast<int>(0xffffffff != idIndex);
	outfile << "FIELD FieldData " << std::to_string(numFields) << "\n";
	// iterate over attributes
	for (int a = 0; a < partioData->numAttributes(); a++)
	{
		if (posIndex == a || idIndex == a)
			continue;

		Partio::ParticleAttribute attr;
		partioData->attributeInfo(a, attr);
		std::string attrNameVTK;
		std::regex_replace(std::back_inserter(attrNameVTK), attr.name.begin(), attr.name.end(), std::regex("\\s+"), "_");
		// write header information
		outfile << attrNameVTK << " " << attr.count << " " << numParticles;
		// write depending on data type
		if (attr.type == Partio::ParticleAttributeType::FLOAT)
		{
			outfile << " float\n";
			// copy from partio data
			std::vector<float> attrData;
			attrData.reserve(partioData->numParticles());
			for (unsigned int i = 0u; i < numParticles; i++)
				attrData.emplace_back(*partioData->data<float>(attr, i));
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				swapByteOrder(&attrData[i]);
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(float));
		}
		else if (attr.type == Partio::ParticleAttributeType::VECTOR)
		{
			outfile << " float\n";
			// copy from partio data
			std::vector<Vector3f> attrData;
			attrData.reserve(partioData->numParticles());
			for (unsigned int i = 0u; i < numParticles; i++)
				attrData.emplace_back(partioData->data<float>(attr, i));
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				for (unsigned int c = 0; c < 3; c++)
					swapByteOrder(&attrData[i][c]);
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData[0].data()), 3 * numParticles * sizeof(float));
		}
		else if (attr.type == Partio::ParticleAttributeType::INT)
		{
			outfile << " int\n";
			// copy from partio data
			std::vector<int> attrData;
			attrData.reserve(partioData->numParticles());
			for (unsigned int i = 0u; i < numParticles; i++)
				attrData.emplace_back(*partioData->data<int>(attr, i));
			// swap endianess
			for (unsigned int i = 0; i < numParticles; i++)
				swapByteOrder(&attrData[i]);
			// export to vtk
			outfile.write(reinterpret_cast<char*>(attrData.data()), numParticles * sizeof(int));
		}
		else
		{
			LOG_WARN << "Skipping attribute " << attr.name << ", because it is of unsupported type " << (attr.type == Partio::ParticleAttributeType::INDEXEDSTR ? "INDEXEDSTR" : "NONE") << "\n";
			continue;
		}
		// end of block
		outfile << "\n";
	}
	outfile.close();
	STOP_TIMING_AVG;
}
