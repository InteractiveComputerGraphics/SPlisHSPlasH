#ifndef __Elasticity_Peer2018_h__
#define __Elasticity_Peer2018_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ElasticityBase.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"

namespace SPH
{
	/** \brief This class implements the implicit SPH formulation for 
	* incompressible linearly elastic solids introduced
	* by Peer et al. [PGBT17].
	*
	* References:
	* - [PGBT17] Andreas Peer, Christoph Gissler, Stefan Band, and Matthias Teschner. An implicit SPH formulation for incompressible linearly elastic solids. Computer Graphics Forum, 2017. URL: http://dx.doi.org/10.1111/cgf.13317
	*/
	class Elasticity_Peer2018 : public ElasticityBase
	{
	protected:
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> Solver;

		// initial particle indices, used to access their original positions
		std::vector<unsigned int> m_current_to_initial_index;
		std::vector<unsigned int> m_initial_to_current_index;
		// initial particle neighborhood
		std::vector<std::vector<unsigned int>> m_initialNeighbors;
		// volumes in rest configuration
		std::vector<Real> m_restVolumes;
		std::vector<Matrix3r> m_rotations;
		std::vector<Vector6r> m_stress;
		std::vector<Matrix3r> m_L;
		std::vector<Matrix3r> m_RL;
		std::vector<Matrix3r> m_F;
		unsigned int m_iterations;
		unsigned int m_maxIter;
		Real m_maxError;
		Real m_alpha;
		Solver m_solver;

		void initValues();
		void computeMatrixL();
		void computeRotations();
		void computeRHS(VectorXr & rhs);	

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


	public:
		static int ITERATIONS;
		static int MAX_ITERATIONS;
		static int MAX_ERROR;
		static int ALPHA;

		Elasticity_Peer2018(FluidModel *model);
		virtual ~Elasticity_Peer2018(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Elasticity_Peer2018(model); }

		virtual void step();
		virtual void reset();
		virtual void performNeighborhoodSearchSort();

		virtual void saveState(BinaryFileWriter &binWriter);
		virtual void loadState(BinaryFileReader &binReader);

		static void matrixVecProd(const Real* vec, Real *result, void *userData);
	};
}

#endif
