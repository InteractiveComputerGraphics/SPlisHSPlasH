#ifndef __BinaryFileWriter_H__
#define __BinaryFileWriter_H__

#include <iostream>
#include <fstream>

#include "SPlisHSPlasH/Common.h"

namespace SPH
{
	class BinaryFileWriter
	{
	public:
		static std::ofstream m_file;

	public:
		static bool openMeshFile(const char *fileName);
		static void closeMeshFile();
		static void writeBuffer(const char *buffer, size_t size);
		static void writeMatrix3f(const Eigen::Matrix3f &m);
		static void writeMatrix3d(const Eigen::Matrix3d &m);
		static void writeVector3f(const Eigen::Vector3f &v);
		static void writeVector3d(const Eigen::Vector3d &v);
		static void writeFloat(const float v);
		static void writeDouble(const double v);
		static void writeString(const std::string &str);
		static void writeBool(const bool b);
		static void writeInt(const int i);
		static void writeUInt(const unsigned int i);
		static void writeChar(const char c);	
	};
}

#endif