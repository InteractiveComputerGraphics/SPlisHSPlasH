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
	template <typename _Scalar>
	class JacobiPreconditioner
	{
	public:
		typedef typename SystemMatrixType::StorageIndex StorageIndex;
		typedef void(*DiagonalMatrixElementFct) (const unsigned int, Real&, void *);

		enum {
			ColsAtCompileTime = Eigen::Dynamic,
			MaxColsAtCompileTime = Eigen::Dynamic
		};

		JacobiPreconditioner() {}

		void init(const unsigned int dim, DiagonalMatrixElementFct fct, void *userData)
		{
			m_dim = dim; m_diagonalElementFct = fct; m_userData = userData;
		}

		Eigen::Index rows() const { return m_dim; }
		Eigen::Index cols() const { return m_dim; }

		Eigen::ComputationInfo info() { return Eigen::Success; }

		template<typename MatType>
		JacobiPreconditioner& analyzePattern(const MatType&) { return *this; }

		template<typename MatType>
		JacobiPreconditioner& factorize(const MatType& mat) { return *this; }

		template<typename MatType>
		JacobiPreconditioner& compute(const MatType& mat) 
		{ 
			m_invDiag.resize(m_dim);
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static) 
				for (int i = 0; i < (int)m_dim; i++)
				{
					Real res;
					m_diagonalElementFct(i, res, m_userData);
					m_invDiag[i] = 1.0 / res;
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
		inline const Eigen::Solve<JacobiPreconditioner, Rhs> solve(const Eigen::MatrixBase<Rhs>& b) const
		{
			return Eigen::Solve<JacobiPreconditioner, Rhs>(*this, b.derived());
		}

	protected:
		unsigned int m_dim;
		/** diagonal matrix element callback */
		DiagonalMatrixElementFct m_diagonalElementFct;
		void *m_userData;
		Eigen::VectorXd m_invDiag;
	};
}

namespace Eigen
{
	namespace internal
	{
		using namespace SPH;

		/** Implementation of the matrix-free matrix vector product  */
		template<typename Rhs>
		struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for generic matrix-vector
			: generic_product_impl_base<MatrixReplacement, Rhs, generic_product_impl<MatrixReplacement, Rhs> >
		{
			typedef typename Product<MatrixReplacement, Rhs>::Scalar Scalar;

			template<typename Dest>
			static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
			{
				// This method should implement "dst += alpha * lhs * rhs" inplace,
				// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
				assert(alpha == Scalar(1) && "scaling is not implemented");

				const Real *vec = &rhs(0);
				Real *res = &dst(0);
				MatrixReplacement& lhs_ = const_cast<MatrixReplacement&>(lhs);
				lhs_.getMatrixVecProdFct()(vec, res, lhs_.getUserData());
			}
		};
	}
}


#endif
