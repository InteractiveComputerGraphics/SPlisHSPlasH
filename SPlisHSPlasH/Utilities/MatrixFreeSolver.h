#ifndef __MatrixFreeSolver_H__
#define __MatrixFreeSolver_H__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using SystemMatrixType = Eigen::SparseMatrix<Real>;

namespace SPH
{
	class MatrixReplacement;
}

namespace Eigen
{
	namespace internal
	{
		template<> struct traits<SPH::MatrixReplacement> : public Eigen::internal::traits<SystemMatrixType> {};
	}
}

namespace SPH
{
	/** Replacement of the matrix in the linear system which is required for a
	* matrix-free solver. */
	class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement>
	{
	public:
		// Required typedefs, constants, and method:
		typedef Real Scalar;
		typedef Real RealScalar;
		typedef int StorageIndex;
		typedef void(*MatrixVecProdFct) (const Real*, Real*, void *);

		enum
		{
			ColsAtCompileTime = Eigen::Dynamic,
			MaxColsAtCompileTime = Eigen::Dynamic,
			IsRowMajor = false
		};

		Index rows() const { return m_dim; }
		Index cols() const { return m_dim; }

		template<typename Rhs>
		Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const
		{
			return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
		}

		MatrixReplacement(const unsigned int dim, MatrixVecProdFct fct, void *userData) : m_dim(dim), m_matrixVecProdFct(fct), m_userData(userData) {}
		void * getUserData() { return m_userData; }
		MatrixVecProdFct getMatrixVecProdFct() { return m_matrixVecProdFct; }

	protected:
		unsigned int m_dim;
		void *m_userData;
		/** matrix vector product callback */
		MatrixVecProdFct m_matrixVecProdFct;
	};


	/** Matrix-free Jacobi preconditioner  */
	class JacobiPreconditioner1D
	{
	public:
		typedef typename SystemMatrixType::StorageIndex StorageIndex;
		typedef void(*DiagonalMatrixElementFct) (const unsigned int, Real&, void *);

		enum {
			ColsAtCompileTime = Eigen::Dynamic,
			MaxColsAtCompileTime = Eigen::Dynamic
		};

		JacobiPreconditioner1D() {}

		void init(const unsigned int dim, DiagonalMatrixElementFct fct, void *userData)
		{
			m_dim = dim; m_diagonalElementFct = fct; m_userData = userData;
		}

		Eigen::Index rows() const { return m_dim; }
		Eigen::Index cols() const { return m_dim; }

		Eigen::ComputationInfo info() { return Eigen::Success; }

		template<typename MatType>
		JacobiPreconditioner1D& analyzePattern(const MatType&) { return *this; }

		template<typename MatType>
		JacobiPreconditioner1D& factorize(const MatType& mat) { return *this; }

		template<typename MatType>
		JacobiPreconditioner1D& compute(const MatType& mat) 
		{ 
			m_invDiag.resize(m_dim);
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static) 
				for (int i = 0; i < (int)m_dim; i++)
				{
					Real res;
					m_diagonalElementFct(i, res, m_userData);
					m_invDiag[i] = static_cast<Real>(1.0) / res;
				}
			}
			return *this; 
		}

		template<typename Rhs, typename Dest>
		void _solve_impl(const Rhs& b, Dest& x) const
		{
			x = m_invDiag.array() * b.array();
		}

		template<typename Rhs>
		inline const Eigen::Solve<JacobiPreconditioner1D, Rhs> solve(const Eigen::MatrixBase<Rhs>& b) const
		{
			return Eigen::Solve<JacobiPreconditioner1D, Rhs>(*this, b.derived());
		}

	protected:
		unsigned int m_dim;
		/** diagonal matrix element callback */
		DiagonalMatrixElementFct m_diagonalElementFct;
		void *m_userData;
		VectorXr m_invDiag;
	};

	/** Matrix-free Jacobi preconditioner  */
	class JacobiPreconditioner3D
	{
	public:
		typedef typename SystemMatrixType::StorageIndex StorageIndex;
		typedef void(*DiagonalMatrixElementFct) (const unsigned int, Vector3r&, void *);
		
		enum {
			ColsAtCompileTime = Eigen::Dynamic,
			MaxColsAtCompileTime = Eigen::Dynamic
		};

		JacobiPreconditioner3D() {}

		void init(const unsigned int dim, DiagonalMatrixElementFct fct, void *userData)
		{
			m_dim = dim; m_diagonalElementFct = fct; m_userData = userData;
		}

		Eigen::Index rows() const { return 3*m_dim; }
		Eigen::Index cols() const { return 3*m_dim; }

		Eigen::ComputationInfo info() { return Eigen::Success; }

		template<typename MatType>
		JacobiPreconditioner3D& analyzePattern(const MatType&) { return *this; }

		template<typename MatType>
		JacobiPreconditioner3D& factorize(const MatType& mat) { return *this; }

		template<typename MatType>
		JacobiPreconditioner3D& compute(const MatType& mat)
		{
			m_invDiag.resize(m_dim*3);
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static) 
				for (int i = 0; i < (int)m_dim; i++)
				{
					Vector3r res;
					m_diagonalElementFct(i, res, m_userData);
					m_invDiag[3*i] = static_cast<Real>(1.0) / res[0];
					m_invDiag[3*i+1] = static_cast<Real>(1.0) / res[1];
					m_invDiag[3*i+2] = static_cast<Real>(1.0) / res[2];
				}
			}
			return *this;
		}

		template<typename Rhs, typename Dest>
		void _solve_impl(const Rhs& b, Dest& x) const
		{
			x = m_invDiag.array() * b.array();
		}

		template<typename Rhs>
		inline const Eigen::Solve<JacobiPreconditioner3D, Rhs> solve(const Eigen::MatrixBase<Rhs>& b) const
		{
			return Eigen::Solve<JacobiPreconditioner3D, Rhs>(*this, b.derived());
		}

	protected:
		unsigned int m_dim;
		/** diagonal matrix element callback */
		DiagonalMatrixElementFct m_diagonalElementFct;
		void *m_userData;
		VectorXr m_invDiag;
	};

	/** Matrix-free 3x3 block Jacobi preconditioner  */
	class BlockJacobiPreconditioner3D
	{
	public:
		typedef typename SystemMatrixType::StorageIndex StorageIndex;
		typedef void(*DiagonalMatrixElementFct) (const unsigned int, Matrix3r&, void *);

		enum {
			ColsAtCompileTime = Eigen::Dynamic,
			MaxColsAtCompileTime = Eigen::Dynamic
		};

		BlockJacobiPreconditioner3D() {}

		void init(const unsigned int dim, DiagonalMatrixElementFct fct, void *userData)
		{
			m_dim = dim; m_diagonalElementFct = fct; m_userData = userData;
		}

		Eigen::Index rows() const { return 3 * m_dim; }
		Eigen::Index cols() const { return 3 * m_dim; }

		Eigen::ComputationInfo info() { return Eigen::Success; }

		template<typename MatType>
		BlockJacobiPreconditioner3D& analyzePattern(const MatType&) { return *this; }

		template<typename MatType>
		BlockJacobiPreconditioner3D& factorize(const MatType& mat) { return *this; }

		template<typename MatType>
		BlockJacobiPreconditioner3D& compute(const MatType& mat)
		{
			m_invDiag.resize(m_dim);
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static) 
				for (int i = 0; i < (int)m_dim; i++)
				{
					Matrix3r res;
					m_diagonalElementFct(i, res, m_userData);
					m_invDiag[i] = res.inverse();
				}
			}
			return *this;
		}

		template<typename Rhs, typename Dest>
		void _solve_impl(const Rhs& b, Dest& x) const
		{
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static) 
				for (int i = 0; i < (int)m_dim; i++)
				{
					static_cast<VectorXr&>(x).block<3, 1>(3 * i, 0) = m_invDiag[i] * static_cast<const VectorXr&>(b).block<3, 1>(3 * i, 0);
				}
			}
		}

		template<typename Rhs>
		inline const Eigen::Solve<BlockJacobiPreconditioner3D, Rhs> solve(const Eigen::MatrixBase<Rhs>& b) const
		{
			return Eigen::Solve<BlockJacobiPreconditioner3D, Rhs>(*this, b.derived());
		}

	protected:
		unsigned int m_dim;
		/** diagonal matrix element callback */
		DiagonalMatrixElementFct m_diagonalElementFct;
		void *m_userData;
		std::vector<Matrix3r> m_invDiag;
	};
}

namespace Eigen
{
	namespace internal
	{
		/** Implementation of the matrix-free matrix vector product  */
		template<typename Rhs>
		struct generic_product_impl<SPH::MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for generic matrix-vector
			: generic_product_impl_base<SPH::MatrixReplacement, Rhs, generic_product_impl<SPH::MatrixReplacement, Rhs> >
		{
			typedef typename Product<SPH::MatrixReplacement, Rhs>::Scalar Scalar;

			template<typename Dest>
			static void scaleAndAddTo(Dest& dst, const SPH::MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
			{
				// This method should implement "dst += alpha * lhs * rhs" inplace,
				// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
				assert(alpha == Scalar(1) && "scaling is not implemented");

				const Real *vec = &rhs(0);
				Real *res = &dst(0);
				SPH::MatrixReplacement& lhs_ = const_cast<SPH::MatrixReplacement&>(lhs);
				lhs_.getMatrixVecProdFct()(vec, res, lhs_.getUserData());
			}
		};
	}
}


#endif
