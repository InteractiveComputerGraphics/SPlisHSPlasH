#ifndef __Viscosity_Peer2016_h__
#define __Viscosity_Peer2016_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"


namespace SPH
{
	/** \brief This class implements the implicit simulation method for
	* viscous fluids introduced
	* by Peer and Teschner [PT16].
	*
	* References:
	* - [PT16] Andreas Peer and Matthias Teschner. Prescribed Velocity Gradients for Highly Viscous SPH Fluids with Vorticity Diffusion. IEEE Transactions on Visualization and Computer Graphics, 2016. URL: https://doi.org/10.1109/TVCG.2016.2636144
	*/
	class Viscosity_Peer2016 : public ViscosityBase
	{
	protected: 
		std::vector<Real> m_density;
		std::vector<Matrix3r> m_targetNablaV;
		std::vector<Vector3r> m_omega;
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner1D> Solver;
		Solver m_solverV;
		Solver m_solverOmega;
		unsigned int m_iterationsV;
		unsigned int m_iterationsOmega;
		unsigned int m_maxIterV;
		Real m_maxErrorV;
		unsigned int m_maxIterOmega;
		Real m_maxErrorOmega;

		virtual void initParameters();
		void computeDensities();

	public:
		static int ITERATIONS_V;
		static int ITERATIONS_OMEGA;
		static int MAX_ITERATIONS_V;
		static int MAX_ERROR_V;
		static int MAX_ITERATIONS_OMEGA;
		static int MAX_ERROR_OMEGA;

		Viscosity_Peer2016(FluidModel *model);
		virtual ~Viscosity_Peer2016(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Viscosity_Peer2016(model); }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		static void matrixVecProdV(const Real* vec, Real *result, void *userData);
		FORCE_INLINE static void diagonalMatrixElementV(const unsigned int row, Real &result, void *userData);

		static void matrixVecProdOmega(const Real* vec, Real *result, void *userData);
		FORCE_INLINE static void diagonalMatrixElementOmega(const unsigned int row, Real &result, void *userData);

		FORCE_INLINE const Matrix3r& getTargetNablaV(const unsigned int i) const
		{
			return m_targetNablaV[i];
		}

		FORCE_INLINE Matrix3r& getTargetNablaV(const unsigned int i)
		{
			return m_targetNablaV[i];
		}

		FORCE_INLINE void setTargetNablaV(const unsigned int i, const Matrix3r &val)
		{
			m_targetNablaV[i] = val;
		}

		FORCE_INLINE const Vector3r& getOmega(const unsigned int i) const
		{
			return m_omega[i];
		}

		FORCE_INLINE Vector3r& getOmega(const unsigned int i)
		{
			return m_omega[i];
		}

		FORCE_INLINE void setOmega(const unsigned int i, const Vector3r &val)
		{
			m_omega[i] = val;
		}
	};
}

#endif
