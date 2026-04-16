#ifndef __Elasticity_Kee2023_h__
#define __Elasticity_Kee2023_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"
#if USE_AVX
#include "SPlisHSPlasH/Utilities/AVX_math.h"
#include "SPlisHSPlasH/Utilities/CholeskyAVXSolver.h"
#endif


namespace SPH
{
	/** \brief This class implements the elasticity solver
	* by Kee et al. [K23].
	*
	* References:
	* - [K23] Kee et al..
	* Kee 2023.
	*
	*/
	class Elasticity_Kee2023 : public NonPressureForceBase
	{
	protected:

		struct Factorization
		{
			Real m_dt;
			Real m_mu;
			Eigen::SparseMatrix<Real, Eigen::RowMajor> m_DT_K;
			Eigen::SparseMatrix<Real, Eigen::RowMajor> m_D;
			Eigen::SparseMatrix<Real, Eigen::ColMajor> m_matHTH;

#ifdef USE_AVX
			CholeskyAVXSolver *m_cholesky;
			Factorization() { m_cholesky = nullptr; }
			~Factorization() { delete m_cholesky; }
#else
			Factorization() {}
			~Factorization() {}
			// L-BFGS prefactored Cholesky (N×N, constant proxy)
			Eigen::SparseMatrix<Real, Eigen::ColMajor> m_matL;
			Eigen::SparseMatrix<Real, Eigen::ColMajor> m_matLT;
			Eigen::VectorXi m_permInd;
			Eigen::VectorXi m_permInvInd;
#endif
		};

		struct ElasticObject
		{
			std::string m_md5;
			std::vector<unsigned int> m_particleIndices;
			unsigned int m_nFixed;
			std::shared_ptr<Factorization> m_factorization;

			// Newton: per-particle 9×9 Hessian and block-diagonal preconditioner.
			// Same type in both builds (SVD is always scalar).
			std::vector<Eigen::Matrix<Real, 9, 9>> m_hessian9x9;
			std::vector<Matrix3r, Eigen::aligned_allocator<Matrix3r>> m_pcg_precond;

#ifdef USE_AVX
			// AVX state (Scalarf8, 3 active lanes x/y/z).
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_f_avx;          // F = D·xk workspace (3*numParticles)
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_xk_avx;         // current iterate
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_xTilde_avx;     // inertial target
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_dx_avx;         // solver step (also LLT solve RHS/result)
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_gradient_avx;   // ∇E at xk

			// Newton CG workspace (Scalarf8, coord-packed x/y/z in lanes 0/1/2)
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_pcg_r_avx;
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_pcg_p_avx;
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_pcg_Ap_avx;
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_pcg_z_avx;

			// L-BFGS secant history
			std::vector<std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>>> m_lbfgs_s_avx;
			std::vector<std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>>> m_lbfgs_y_avx;
			std::vector<Real> m_lbfgs_rho;
			std::vector<Real> m_lbfgs_alpha;
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_lbfgs_last_sol_avx;
			std::vector<Scalarf8, AlignmentAllocator<Scalarf8, 32>> m_lbfgs_q_avx;
			int m_lbfgs_count = 0;
#else
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_f;        // F = D·xk workspace
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_xk;       // current iterate
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_xTilde;   // inertial target
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_dx;       // Newton/LBFGS step (also LLT solve RHS/result)
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_gradient; // ∇E at xk

			// Newton PCG workspace
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_pcg_r;
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_pcg_p;
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_pcg_Ap;
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_pcg_z;

			// L-BFGS secant history
			std::vector<std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>>> m_lbfgs_s;
			std::vector<std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>>> m_lbfgs_y;
			std::vector<Real> m_lbfgs_rho;
			std::vector<Real> m_lbfgs_alpha;
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_lbfgs_last_sol;
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_lbfgs_q;
			int m_lbfgs_count = 0;

			// Permutation workspace for scalar LLT (manual forward/backward sub)
			std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_dx_perm;
#endif

			ElasticObject()  { m_factorization = nullptr; }
			~ElasticObject() { m_factorization = nullptr; }
		};

		Real m_youngsModulus;
		Real m_poissonRatio;
		Vector3r m_fixedBoxMin;
		Vector3r m_fixedBoxMax;
		Vector3r m_fixedBox2Min;
		Vector3r m_fixedBox2Max;

		// initial particle indices, used to access their original positions
		std::vector<unsigned int> m_current_to_initial_index;
		std::vector<unsigned int> m_initial_to_current_index;
		// initial particle neighborhood
		std::vector<std::vector<unsigned int>> m_initialNeighbors;
		// volumes in rest configuration
		std::vector<Real> m_restVolumes;
		std::vector<Matrix3r> m_rotations;
		std::vector<Real> m_stress;
		std::vector<int> m_fixedGroupId;		// 0: free, 1: box1, 2: box2
		std::vector<Matrix3r> m_L;
		std::vector<Matrix3r> m_F;
		std::vector<Matrix3r> m_PL;
		Real m_alpha;
		int m_maxNeighbors;
		int m_solverType;			// 0: Newton, 1: LBFGS
		int m_lbfgsWindowSize;
		int m_materialType;			// 0: Stable Neo-Hookean, 1: Co-rotated
		int m_maxIter;
		Real m_maxError;
		int m_maxIterCG;      // max CG iterations for Newton linear solve
		Real m_tolCG;         // CG convergence tolerance
		int m_maxLSIter;
		Real m_lsArmijoParam;
		Real m_lsBeta;
		bool m_useLineSearch;
		unsigned int m_totalNeighbors;
		std::vector<ElasticObject*> m_objects;
		Real m_lambda;
		Real m_mu;

		// Precomputed V_j * gradW(xi0 - xj0) per neighbor pair.
		// Scalar format (always available, used by Newton CG).
		std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_precomp_V_gradW;
		std::vector<unsigned int> m_precomputed_indices;
#ifdef USE_AVX
		// AVX-packed format (8 neighbors per entry, used by L-BFGS force loops).
		std::vector<Vector3f8, Eigen::aligned_allocator<Vector3f8>> m_precomp_V_gradW8;
		std::vector<unsigned int> m_precomputed_indices8;
#endif

#ifdef USE_AVX
		typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> SolverLLT;
#else
		typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> SolverLLT;
#endif

		void determineFixedParticles();
		std::string computeMD5(const unsigned int objIndex);
		void initValues();
		void initSystem();
		void initFactorization(std::shared_ptr<Factorization> factorization, std::vector<unsigned int> &particleIndices, const unsigned int nFixed, const Real dt, const Real mu);
		void findObjects();
		void computeMatrixL();
		void precomputeValues();

		void stepElasticitySolver();

		void computeXTilde(ElasticObject* obj);
		void updateVelocity(ElasticObject* obj, Real fdt);
		Real computeEnergy(ElasticObject* obj);
		Real computePsi(const Matrix3r& F, const Matrix3r& R) const;
		Real computeEnergyAndGradient(ElasticObject* obj);
		void computeHessian(ElasticObject* obj);
		void computeCorotatedHessian9x9(ElasticObject* obj);
		void computeStableNeoHookeanHessian9x9(ElasticObject* obj);
		void computeNewtonPreconditioner(ElasticObject* obj);
		void newtonMatvec(ElasticObject* obj);
		int matFreePCG(ElasticObject* obj);
		Real newtonSolve(ElasticObject* obj, int& cgIter);
		void prefactorizedLLTSolve(ElasticObject* obj);
		Real lbfgsSolve(ElasticObject* obj);
		Real lineSearch(ElasticObject* obj, Real energy, int& lsIter);

		Matrix3r computeP(const Matrix3r& F, const Matrix3r& R) const;

		virtual void initParameters();
		/** This function is called after the simulation scene is loaded and all
		* parameters are initialized. While reading a scene file several parameters
		* can change. The deferred init function should initialize all values which
		* depend on these parameters.
		*/
		virtual void deferredInit();

		//////////////////////////////////////////////////////////////////////////
		// multiplication of symmetric matrix, represented by a 6D vector, and a
		// 3D vector
		//////////////////////////////////////////////////////////////////////////
		FORCE_INLINE void symMatTimesVec(const Vector6r & M, const Vector3r & v, Vector3r &res)
		{
			res[0] = M[0] * v[0] + M[3] * v[1] + M[4] * v[2];
			res[1] = M[3] * v[0] + M[1] * v[1] + M[5] * v[2];
			res[2] = M[4] * v[0] + M[5] * v[1] + M[2] * v[2];
		}

		void rotationMatricesToAVXQuaternions();

	public:
		static std::string METHOD_NAME;
		static int YOUNGS_MODULUS;
		static int POISSON_RATIO;
		static int FIXED_BOX_MIN;
		static int FIXED_BOX_MAX;
		static int FIXED_BOX2_MIN;
		static int FIXED_BOX2_MAX;
		static int ALPHA;
		static int MAX_NEIGHBORS;
		static int SOLVER_TYPE;
		static int LBFGS_WINDOW_SIZE;
		static int MATERIAL_TYPE;
		static int MAX_ITER;
		static int MAX_ERROR;
		static int MAX_ITER_CG;
		static int TOL_CG;
		static int MAX_LS_ITER;
		static int LS_ARMIJO_PARAM;
		static int LS_BETA;
		static int USE_LINE_SEARCH;

		static int ENUM_SOLVER_NEWTON;
		static int ENUM_SOLVER_LBFGS;
		static int ENUM_MATERIAL_STABLE_NEOHOOKEAN;
		static int ENUM_MATERIAL_COROTATED;

		Elasticity_Kee2023(FluidModel *model);
		virtual ~Elasticity_Kee2023(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Elasticity_Kee2023(model); }
		virtual std::string getMethodName() { return METHOD_NAME; }

		virtual void step();
		virtual void reset();
		virtual void performNeighborhoodSearchSort();

		virtual void saveState(BinaryFileWriter &binWriter);
		virtual void loadState(BinaryFileReader &binReader);
	};
}

#endif
