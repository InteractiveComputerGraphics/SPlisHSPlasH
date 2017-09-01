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
	* by Peer and Teschner \cite Peer2016.
	*/
	class Viscosity_Peer2016 : public ViscosityBase
	{
	protected: 
		std::vector<Matrix3r> m_targetNablaV;
		std::vector<Vector3r> m_omega;
		unsigned int m_maxIter;
		Real m_maxError;
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner<Real>> Solver;
		Solver m_solverV;
		Solver m_solverOmega;

	public:
		Viscosity_Peer2016(FluidModel *model);
		virtual ~Viscosity_Peer2016(void);

		unsigned int getMaxIter() const { return m_maxIter; }
		void setMaxIter(unsigned int val) { m_maxIter = val; }
		Real getMaxError() const { return m_maxError; }
		void setMaxError(Real val) { m_maxError = val; }

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
