#include "CholeskyAVXSolver.h"
#include "Utilities/Logger.h"

using namespace SPH;

CholeskyAVXSolver::CholeskyAVXSolver(const Eigen::SparseMatrix<float>& lhs)
{
	this->_init(lhs.cast<double>());
}

CholeskyAVXSolver::CholeskyAVXSolver(const Eigen::SparseMatrix<double>& lhs)
{
	this->_init(lhs);
}

void CholeskyAVXSolver::solve(float* solution, const float* rhs, const int stride)
{
	// Common
	const CholeskySparseMatrixfAVX& L = this->L;
	const CholeskySparseMatrixfAVX& LT = this->LT;

	// AVX rhs
	std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> rhs_avx(this->n_rhs_lines);  // mem request for thread safe
	rhs_avx.back() = Scalarf8(0.0f); // For potential padding

	union AVXfloatUnion {
		Scalarf8* lines;
		float* f;
	};
	AVXfloatUnion avx_float_union = { rhs_avx.data() };
	float* rhs_scalar = avx_float_union.f;

	// Copy and permute input to the extended rhs
	for (int i = 0; i < this->ndofs; i++) {
		rhs_scalar[this->perm[i]] = rhs[stride*i];
	}

	// Solve linear system
	//// Forward substitution
	for (int row = 0; row < this->ndofs; row++) {
		Scalarf8 row_sum_avx = Scalarf8(0.0f);
		for (int col_offset = L.rows_offset[row]; col_offset < L.rows_offset[row + 1]; col_offset++) {
			const int& col = L.cols[col_offset];
			const Scalarf8& val = L.vals[col_offset];
			row_sum_avx += val * rhs_avx[col];
		}

		// Sum items in line
		//float row_sum = row_sum_avx.hadd();
		const float row_sum = row_sum_avx.reduce();

		// Substitute
		rhs_scalar[row] = (rhs_scalar[row] - row_sum) * this->L.diagonal_inv[row];  // L and LT have the same diagonal
	}

	//// Backward substitution
	for (int row = this->ndofs - 1; row >= 0; row--) {
		Scalarf8 row_sum_avx = Scalarf8(0.0f);
		for (int col_offset = LT.rows_offset[row]; col_offset < LT.rows_offset[row + 1]; col_offset++) {
			const int& col = LT.cols[col_offset];
			const Scalarf8& val = LT.vals[col_offset];
			row_sum_avx += val * rhs_avx[col];
		}

		// Sum items in line
		//float row_sum = row_sum_avx.hadd();
		const float row_sum = row_sum_avx.reduce();

		// Substitute
		rhs_scalar[row] = (rhs_scalar[row] - row_sum) * this->LT.diagonal_inv[row];  // L and LT have the same diagonal
	}

	// Copy and permute input to the extended rhs
	for (int i = 0; i < this->ndofs; i++) {
		solution[stride*this->perm_inv[i]] = rhs_scalar[i];
	}
}

void CholeskyAVXSolver::_init(const Eigen::SparseMatrix<double>& lhs)
{
	// Logic
	this->ndofs = (int)lhs.rows();
	const int extra_dofs = ndofs % 8;
	this->n_rhs_lines = (extra_dofs > 0) ? ndofs / 8 + 1 : ndofs / 8;

	// TODO: Once this is working do not store this tmp variables
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> llt;
	llt.compute(lhs);

	if (llt.info() != Eigen::Success) {
		std::cout << "Cholesky decomposition failed. Error: " << llt.info();
		exit(-1);
	}

	Eigen::SparseMatrix<float, Eigen::RowMajor> L_eigen = Eigen::SparseMatrix<float, Eigen::RowMajor>(llt.matrixL().cast<float>());
	Eigen::SparseMatrix<float, Eigen::RowMajor> LT_eigen = Eigen::SparseMatrix<float, Eigen::RowMajor>(llt.matrixU().cast<float>());

	LOG_INFO << "Non zero elements (L): " << L_eigen.nonZeros();

	auto& eigen_permutation = llt.permutationP().indices();
	this->perm.insert(this->perm.end(), &eigen_permutation[0], &eigen_permutation[0] + this->ndofs);

	auto& eigen_permutation_inv = llt.permutationPinv().indices();
	this->perm_inv.insert(this->perm_inv.end(), &eigen_permutation_inv[0], &eigen_permutation_inv[0] + this->ndofs);

	this->L = CholeskySparseMatrixfAVX(L_eigen);
	this->LT = CholeskySparseMatrixfAVX(LT_eigen);
}

CholeskySparseMatrixfAVX::CholeskySparseMatrixfAVX(const Eigen::SparseMatrix<float, Eigen::RowMajor>& lhs)
{
	const int ndofs = (int)lhs.rows();
	const int extra_dofs = ndofs % 8;
	this->n_rhs_lines = (extra_dofs > 0) ? ndofs/8 + 1 : ndofs/8;
	const int ndofs_ceil = 8 * this->n_rhs_lines;

	this->diagonal_inv.resize(ndofs);
	this->rows_offset.push_back(0);
	std::vector<float> dense_row_buffer(ndofs_ceil, 0.0f);
	std::vector<bool> active_lines(this->n_rhs_lines, false);
	for (int row = 0; row < lhs.outerSize(); ++row) {
		// Clear row buffers
		std::fill(dense_row_buffer.begin(), dense_row_buffer.end(), 0.0f);
		std::fill(active_lines.begin(), active_lines.end(), false);
		
		// Fill current row with the data from the sparse matrix
		for (Eigen::SparseMatrix<float, Eigen::RowMajor>::InnerIterator it(lhs, row); it; ++it) {
			if (it.row() == it.col()) {
				this->diagonal_inv[it.row()] = (float)(1.0 / (double)it.value());
			}
			else {
				dense_row_buffer[it.col()] = it.value();
				active_lines[it.col() / 8] = true;
			}
		}

		// Initialize AVX lines
		for (int line_i = 0; line_i < this->n_rhs_lines; line_i++) {
			if (active_lines[line_i]) {
				const int col = 8 * line_i;
				this->cols.push_back(col / 8);
				this->vals.push_back(Scalarf8(dense_row_buffer[col+0], dense_row_buffer[col+1], dense_row_buffer[col+2], dense_row_buffer[col+3], dense_row_buffer[col+4], dense_row_buffer[col+5], dense_row_buffer[col+6], dense_row_buffer[col+7]));
			}
		}
		this->rows_offset.push_back((int)this->vals.size());
	}
}

void CholeskyAVXSolver::save(SPH::BinaryFileWriter& binWriter)
{
	L.save(binWriter);
	LT.save(binWriter);

	binWriter.writeVector(perm);
	binWriter.writeVector(perm_inv);
	binWriter.write(n_rhs_lines);
	binWriter.write(ndofs);
}

void CholeskyAVXSolver::load(SPH::BinaryFileReader& binReader)
{
	L.load(binReader);
	LT.load(binReader);

	binReader.readVector(perm);
	binReader.readVector(perm_inv);
	binReader.read(n_rhs_lines);
	binReader.read(ndofs);
}

void CholeskySparseMatrixfAVX::save(SPH::BinaryFileWriter& binWriter)
{
	binWriter.write(vals.size());
	binWriter.writeBuffer((char*)vals.data(), vals.size() * sizeof(Scalarf8));
	
	binWriter.writeVector(cols);
	binWriter.writeVector(rows_offset);
	binWriter.writeVector(diagonal_inv);
	binWriter.write(n_rhs_lines);
}

void CholeskySparseMatrixfAVX::load(SPH::BinaryFileReader& binReader)
{
	size_t size;
	binReader.read(size);
	vals.resize(size);
	binReader.readBuffer((char*)vals.data(), size * sizeof(Scalarf8));

	binReader.readVector(cols);
	binReader.readVector(rows_offset);
	binReader.readVector(diagonal_inv);
	binReader.read(n_rhs_lines);
}