#ifndef __Viscosity_Peer2015_h__
#define __Viscosity_Peer2015_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"


namespace SPH
{
	/** \brief This class implements the implicit simulation method for
	* viscous fluids introduced
	* by Peer et al. \cite Peer2015.
	*/
	class Viscosity_Peer2015 : public ViscosityBase
	{
	protected: 
		std::vector<Matrix3r> m_targetNablaV;
		unsigned int m_maxIter;
		Real m_maxError;
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner<Real>> Solver;
		Solver m_solver;

	public:
		Viscosity_Peer2015(FluidModel *model);
		virtual ~Viscosity_Peer2015(void);

		unsigned int getMaxIter() const { return m_maxIter; }
		void setMaxIter(unsigned int val) { m_maxIter = val; }
		Real getMaxError() const { return m_maxError; }
		void setMaxError(Real val) { m_maxError = val; }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		static void matrixVecProd(const Real* vec, Real *result, void *userData);
		FORCE_INLINE static void diagonalMatrixElement(const unsigned int row, Real &result, void *userData);

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
	};
}

#endif
