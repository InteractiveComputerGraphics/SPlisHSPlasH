#ifndef __TimeStepPF_h__
#define __TimeStepPF_h__

#include "SimulationDataPF.h"
#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"

// since all diagonal blocks are 3x3 diagonal matrices, a diagonal preconditioner does suffice
#define PD_USE_DIAGONAL_PRECONDITIONER

namespace SPH
{
	/** \brief This class implements the Projective Fluids approach introduced
	* by Weiler, Koschier and Bender [WKB16].
	*
	* References:
	* - [WKB16] Marcel Weiler, Dan Koschier, and Jan Bender. Projective fluids. In Proceedings of the 9th International Conference on Motion in Games, MIG '16, 79-84. New York, NY, USA, 2016. ACM. URL: http://doi.acm.org/10.1145/2994258.2994282
	*/
	class TimeStepPF : public TimeStep
	{
	protected:
		using VectorXr = Eigen::Matrix<Real, -1, 1>;
		using VectorXrMap = Eigen::Map<VectorXr>;

#ifdef PD_USE_DIAGONAL_PRECONDITIONER
		using Solver = Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner3D>;
		FORCE_INLINE static void diagonalMatrixElement(const unsigned int row, Vector3r &result, void *userData);
		void preparePreconditioner();
#else
		using Solver = Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner>;
#endif

		SimulationDataPF m_simulationData;
		Solver m_solver;
		Real m_stiffness;
		unsigned int m_numActiveParticlesTotal;

		void initialGuessForPositions(const unsigned int fluidModelIndex);
		void solvePDConstraints();
		void updatePositionsAndVelocity(const VectorXr & x);
		void addAccellerationToVelocity();

		void matrixFreeRHS(const VectorXr & x, VectorXr & result);

		virtual void performNeighborhoodSearchSort();
		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex) override;

		virtual void initParameters() override;

	public:
		static int STIFFNESS;

		TimeStepPF();
		virtual ~TimeStepPF(void);

		virtual void step()   override;
		virtual void reset()  override;
		virtual void resize() override;

		static void matrixVecProd(const Real* vec, Real *result, void *userData);
	};
}

#endif
