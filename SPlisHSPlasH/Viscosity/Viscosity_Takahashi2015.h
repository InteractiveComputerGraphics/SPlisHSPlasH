#ifndef __Viscosity_Takahashi2015_h__
#define __Viscosity_Takahashi2015_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"


namespace SPH
{
	/** \brief This class implements a variant of the implicit simulation method for
	* viscous fluids introduced by Takahashi et al. \cite Takahashi2015.
	* In the original work of Takahashi et al. the second-ring neighbors are required
	* to create the matrix of the linear system. In contrast we use a meshless 
	* conjugate gradient solver which performs the required matrix-vector multiplication 
	* in two sequential loops. In this way only the one-ring neighbors are required
	* in each loop which increases the performance significantly.\n\n
	* Thanks to Anreas Peer who helped us with the implementation.
	*/
	class Viscosity_Takahashi2015 : public ViscosityBase
	{
	protected: 
		std::vector<Vector3r> m_accel;
		std::vector<Matrix3r> m_viscousStress;
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> Solver;
		Solver m_solver;
		unsigned int m_iterations;
		unsigned int m_maxIter;
		Real m_maxError;

		virtual void initParameters();
		static void computeViscosityAcceleration(Viscosity_Takahashi2015 *visco, const Real* v);

	public:
		static int ITERATIONS;
		static int MAX_ITERATIONS;
		static int MAX_ERROR;

		Viscosity_Takahashi2015(FluidModel *model);
		virtual ~Viscosity_Takahashi2015(void);

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		static void matrixVecProd(const Real* vec, Real *result, void *userData);
		FORCE_INLINE static void diagonalMatrixElement(const unsigned int row, Real &result, void *userData);

		FORCE_INLINE const Matrix3r& getViscousStress(const unsigned int i) const
		{
			return m_viscousStress[i];
		}

		FORCE_INLINE Matrix3r& getViscousStress(const unsigned int i)
		{
			return m_viscousStress[i];
		}

		FORCE_INLINE void setViscousStress(const unsigned int i, const Matrix3r &val)
		{
			m_viscousStress[i] = val;
		}

		FORCE_INLINE const Vector3r& getAccel(const unsigned int i) const
		{
			return m_accel[i];
		}

		FORCE_INLINE Vector3r& getAccel(const unsigned int i)
		{
			return m_accel[i];
		}

		FORCE_INLINE void setAccel(const unsigned int i, const Vector3r &val)
		{
			m_accel[i] = val;
		}
	};
}

#endif
