#ifndef __Viscosity_Weiler2018_h__
#define __Viscosity_Weiler2018_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"

#define USE_BLOCKDIAGONAL_PRECONDITIONER

namespace SPH
{
	/** \brief This class implements the implicit Laplace viscosity method introduced 
	* by Weiler et al. 2018 \cite Weiler2018.
	*/
	class Viscosity_Weiler2018 : public ViscosityBase
	{
	protected:
		Real m_boundaryViscosity;
		unsigned int m_maxIter;
		Real m_maxError;
		unsigned int m_iterations;
		std::vector<Vector3r> m_vDiff;

#ifdef USE_BLOCKDIAGONAL_PRECONDITIONER
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, BlockJacobiPreconditioner3D> Solver;
		FORCE_INLINE static void diagonalMatrixElement(const unsigned int row, Matrix3r &result, void *userData);
#else
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner3D> Solver;
		FORCE_INLINE static void diagonalMatrixElement(const unsigned int row, Vector3r &result, void *userData);
#endif	
		Solver m_solver;


		virtual void initParameters();
		
	public:
		static int ITERATIONS;
		static int MAX_ITERATIONS;
		static int MAX_ERROR;
		static int VISCOSITY_COEFFICIENT_BOUNDARY;

		Viscosity_Weiler2018(FluidModel *model);
		virtual ~Viscosity_Weiler2018(void);

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		static void matrixVecProd(const Real* vec, Real *result, void *userData);
	};
}

#endif
