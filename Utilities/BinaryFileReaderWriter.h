#ifndef __BinaryFileReaderWriter_H__
#define __BinaryFileReaderWriter_H__

#include <iostream>
#include <fstream>

#include "SPlisHSPlasH/Common.h"

namespace SPH
{
	class BinaryFileWriter
	{
	public:
		std::ofstream m_file;

	public:

		bool openFile(const std::string &fileName)
		{
			m_file.open(fileName, std::ios::out | std::ios::binary);
			if (!m_file.is_open())
			{
				std::cout << "Cannot open file.\n";
				return false;
			}
			return true;
		}

		void closeFile()
		{
			m_file.close();
		}

		void writeBuffer(const char *buffer, size_t size)
		{
			m_file.write(buffer, size);
		}

		template<typename T>
		void write(const T &v)
		{
			writeBuffer((char*)&v, sizeof(T));
		}

		void write(const std::string &str)
		{
			write((unsigned int) str.size());
			writeBuffer(str.c_str(), str.size());
		}

		template<typename T>
		void writeMatrix(const T &m)
		{
			writeBuffer((char*)m.data(), m.size() * sizeof(m.data()[0]));
		}
	};
	
	class BinaryFileReader
	{
	public:
		std::ifstream m_file;

	public:

		bool openFile(const std::string &fileName)
		{
			m_file.open(fileName, std::ios::in | std::ios::binary);
			if (!m_file.is_open())
			{
				std::cout << "Cannot open file.\n";
				return false;
			}
			return true;
		}

		void closeFile()
		{
			m_file.close();
		}

		void readBuffer(char *buffer, size_t size)
		{
			m_file.read(buffer, size);
		}

		template<typename T>
		void read(T &v)
		{
			readBuffer((char*)&v, sizeof(T));
		}

		void read(std::string &str)
		{
			unsigned int len;
			read(len);
			char* temp = new char[len + 1];
			readBuffer(temp, len);
			temp[len] = '\0';
			str = temp;
		}

		template<typename T>
		void readMatrix(T &m)
		{
			readBuffer((char*)m.data(), m.size() * sizeof(m.data()[0]));
		}
	};

	
}

#endif