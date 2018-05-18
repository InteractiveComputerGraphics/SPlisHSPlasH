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
	* by Weiler, Koschier and Bender \cite Weiler:2016.
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
		unsigned int m_counter;
		unsigned int m_numActiveParticlesTotal;

		void initialGuessForPositions(const unsigned int fluidModelIndex);
		void solvePDConstraints();
		void updatePositionsAndVelocity(const VectorXr & x);
		void addAccellerationToVelocity();

		void matrixFreeRHS(const VectorXr & x, VectorXr & result);

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();
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
