#ifndef __BinaryFileReaderWriter_H__
#define __BinaryFileReaderWriter_H__

#include <iostream>
#include <fstream>

#include "SPlisHSPlasH/Common.h"
#include <Eigen/Sparse>

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

		template<typename T, int Rows, int Cols>
		void writeMatrixX(const Eigen::Matrix < T, Rows, Cols> & m)
		{
			const Eigen::Index rows = m.rows();
			const Eigen::Index cols = m.cols();
			write(rows);
			write(cols);

			writeBuffer((char*)m.data(), rows * cols* sizeof(T));
		}

		template <typename T, int Options, typename StorageIndex>
		void writeSparseMatrix(Eigen::SparseMatrix<T, Options, StorageIndex>& m) 
		{
			m.makeCompressed();

			const Eigen::Index rows = m.rows();
			const Eigen::Index cols = m.cols();
			const Eigen::Index nnzs = m.nonZeros();
			const Eigen::Index outS = m.outerSize();
			const Eigen::Index innS = m.innerSize();

			write(rows);
			write(cols);
			write(nnzs);
			write(outS);
			write(innS);
				 
			writeBuffer((const char*)(m.valuePtr()), sizeof(T) * nnzs);
			writeBuffer((const char*)(m.outerIndexPtr()), sizeof(StorageIndex) * outS);
			writeBuffer((const char*)(m.innerIndexPtr()), sizeof(StorageIndex) * nnzs);
		}

		template<typename T>
		void writeVector(const std::vector<T>& m)
		{
			write(m.size());
			writeBuffer((char*)m.data(), m.size() * sizeof(T));
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
			char* temp = new char[len + 1u];
			readBuffer(temp, len);
			temp[len] = '\0';
			str = std::string(temp);
			delete[] temp;
		}

		template<typename T>
		void readMatrix(T &m)
		{
			readBuffer((char*)m.data(), m.size() * sizeof(m.data()[0]));
		}

		template<typename T, int Rows, int Cols>
		void readMatrixX(Eigen::Matrix < T, Rows, Cols>& m)
		{
			Eigen::Index rows, cols;
			read(rows);
			read(cols);
			m.resize(rows, cols);

			readBuffer((char*)m.data(), rows * cols * sizeof(T));
		}

		template <typename T, int Options, typename StorageIndex>
		void readSparseMatrix(Eigen::SparseMatrix<T, Options, StorageIndex>& m)
		{
			Eigen::Index rows, cols, nnzs, innS, outS;
			read(rows);
			read(cols);
			read(nnzs);
			read(outS);
			read(innS);

			m.resize(rows, cols);
			m.makeCompressed();
			m.resizeNonZeros(nnzs);

			readBuffer((char*)(m.valuePtr()), sizeof(T) * nnzs);
			readBuffer((char*)(m.outerIndexPtr()), sizeof(StorageIndex) * outS);
			readBuffer((char*)(m.innerIndexPtr()), sizeof(StorageIndex) * nnzs);

			m.finalize();
		}

		template<typename T>
		void readVector(std::vector<T>& m)
		{
			size_t size;
			read(size);
			m.resize(size);
			readBuffer((char*)m.data(), m.size() * sizeof(T));
		}
	};

	
}

#endif