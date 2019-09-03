#include "BinaryFileWriter.h"

using namespace std;
using namespace SPH;

ofstream BinaryFileWriter::m_file;

bool BinaryFileWriter::openMeshFile(const char *fileName)
{
	m_file.open (fileName, ios::out | ios::binary);
	if(!m_file.is_open()) 
	{
		cout << "Cannot open file.";
		return false;
	}
	return true;
}

void BinaryFileWriter::closeMeshFile()
{
	m_file.close();
}

void BinaryFileWriter::writeBuffer(const char *buffer, size_t size)
{
	m_file.write(buffer, size);
}

void BinaryFileWriter::writeBool( const bool b )
{
	writeBuffer((char*) &b, sizeof(bool));
}

void BinaryFileWriter::writeFloat( const float v )
{	
	writeBuffer((char*) &v, sizeof(float));
}

void BinaryFileWriter::writeInt( const int i )
{
	writeBuffer((char*) &i, sizeof(int));
}

void BinaryFileWriter::writeUInt( const unsigned int i )
{
	writeBuffer((char*) &i, sizeof(unsigned int));
}

void BinaryFileWriter::writeChar( const char c )
{
	writeBuffer(&c, sizeof(char));
}

void BinaryFileWriter::writeDouble( const double v )
{
	writeBuffer((char*) &v, sizeof(double));
}

void BinaryFileWriter::writeMatrix3f(const Eigen::Matrix3f &m)
{
	for (unsigned int i=0; i < 3; i++)
	{
		for (unsigned int j=0; j < 3; j++)
		{
			writeFloat(m(i,j));
		}
	}
}

void BinaryFileWriter::writeMatrix3d(const Eigen::Matrix3d &m)
{
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			writeDouble(m(i, j));
		}
	}
}

void BinaryFileWriter::writeVector3f(const Eigen::Vector3f &v)
{
	for (unsigned int j=0; j < 3; j++)
	{
		writeFloat(v[j]);
	}
}

void BinaryFileWriter::writeVector3d(const Eigen::Vector3d &v)
{
	for (unsigned int j = 0; j < 3; j++)
	{
		writeDouble(v[j]);
	}
}


void BinaryFileWriter::writeString( const string &str )
{
	writeBuffer(str.c_str(), str.size());
	writeChar(0);
}

