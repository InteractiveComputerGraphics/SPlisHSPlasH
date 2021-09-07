#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "SPlisHSPlasH/Utilities/AVX_math.h"
#include "Utilities/BinaryFileReaderWriter.h"

namespace SPH
{
	class CholeskySparseMatrixfAVX
	{
	public:
		// Fields
		//// AVX CRS Sparse matrix
		std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> vals;
		std::vector<int> cols;
		std::vector<int> rows_offset;

		//// Other
		std::vector<float> diagonal_inv;
		int n_rhs_lines = -1;

		// Methods
		CholeskySparseMatrixfAVX() {};
		CholeskySparseMatrixfAVX(const Eigen::SparseMatrix<float, Eigen::RowMajor>& lhs);

		void save(SPH::BinaryFileWriter& binWriter);
		void load(SPH::BinaryFileReader& binReader);
	};


	/** Cholesky solver which uses AVX instructions. 
	* Written by José Antonio Fernández-Fernández.
	*/
	class CholeskyAVXSolver
	{
	public:
		// Fields
		CholeskySparseMatrixfAVX L;
		CholeskySparseMatrixfAVX LT;
		std::vector<int> perm;
		std::vector<int> perm_inv;
		int n_rhs_lines = -1;
		int ndofs = -1;

		// Methods
		CholeskyAVXSolver() = default;
		~CholeskyAVXSolver() = default;
		CholeskyAVXSolver(const Eigen::SparseMatrix<float>& lhs);
		CholeskyAVXSolver(const Eigen::SparseMatrix<double>& lhs);

		/*
			Solves the system for one rhs.

			stride indicates the distance between consecutive entries in the rhs. For example,
			to solve a system only for 'y' in a rhs that is stored as 'xyzxyzxyzxyzxyzxyzxyz' 
			and to store the solution in the 'o' spots '_o__o__o__o__o__o__o_', we use:

				CholeskyAVXSolver.solve(solution_begin + 1, rhs_begin + 1, 3);

			The result is written following the same stride.
			IMPORTANT: This method is thread safe.
		*/
		void solve(float* solution, const float* rhs, const int stride);

		void save(SPH::BinaryFileWriter& binWriter);
		void load(SPH::BinaryFileReader& binReader);

	private:
		void _init(const Eigen::SparseMatrix<double>& lhs);
	};
};