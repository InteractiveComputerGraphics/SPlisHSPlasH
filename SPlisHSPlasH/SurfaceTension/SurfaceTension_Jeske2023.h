#ifndef __SurfaceTension_Jeske2023_Surface_Tension_h__
#define __SurfaceTension_Jeske2023_Surface_Tension_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"
#include "Utilities/Logger.h"

namespace SPH
{
	/** \brief This class implements the implicit surface tension method 
	* by Jeske et al. 2023 [JWL+23].
	*
	* References:
	* - [JWL+23] Jeske, Stefan Rhys, Lukas Westhofen, Fabian Löschner, José Antonio Fernández-Fernández, and Jan Bender. “Implicit Surface Tension for SPH Fluid Simulation.” ACM Transactions on Graphics, November 7, 2023. URL: https://doi.org/10.1145/3631936.

	*/
	class SurfaceTension_Jeske2023 : public SurfaceTensionBase
	{
	protected:
        Real m_viscosity;
        Real m_boundaryViscosity;

		unsigned int m_maxIter;
		Real m_maxError;
		unsigned int m_iterations;
		std::vector<Vector3r> m_vDiff;
		std::vector<Real> m_gradRho;
        std::vector<Real> m_surfaceEnergy;
        std::vector<Real> m_color;
        std::vector<Vector3r> m_colorGrad;
        std::vector<Vector3r> m_nonlinearAcc;
        std::vector<Vector3r> m_nonlinearRes;
        std::vector<Vector3r> m_nonlinearGrad;

		Real m_tangentialDistanceFactor;
		bool m_weakPhaseCoupling;
        Real m_xsph;

		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> Solver;

		Solver m_solver;


		virtual void initParameters();
		
	public:
		static int ITERATIONS;
		static int MAX_ITERATIONS;
		static int MAX_ERROR;
        static int VISCOSITY_COEFFICIENT;
		static int VISCOSITY_COEFFICIENT_BOUNDARY;
        static int XSPH;

		SurfaceTension_Jeske2023(FluidModel *model);
		virtual ~SurfaceTension_Jeske2023(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_Jeske2023(model); }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		static void matrixVecProd(const Real* vec, Real *result, void *userData);

		FORCE_INLINE const Vector3r& getVDiff(const unsigned int i) const
		{
			return m_vDiff[i];
		}

		FORCE_INLINE Vector3r& getVDiff(const unsigned int i)
		{
			return m_vDiff[i];
		}

		FORCE_INLINE void setVDiff(const unsigned int i, const Vector3r& val)
		{
			m_vDiff[i] = val;
		}

        FORCE_INLINE const Real& getDensityGrad(const unsigned int i) const
        {
            return m_gradRho[i];
        }

        FORCE_INLINE Real& getDensityGrad(const unsigned int i)
        {
            return m_gradRho[i];
        }

        FORCE_INLINE void setDensityGrad(const unsigned int i, const Real& val)
        {
            m_gradRho[i] = val;
        }

        void computeRHS(VectorXr &b, VectorXr &g);

        void applyForces(const VectorXr &x);

        Real getMaxSolverError(){ return m_maxError; }
        void setMaxSolverError(Real error){ m_maxError = error; }

        bool getWeakCoupling(){ return m_weakPhaseCoupling; }
        void setWeakCoupling(bool val){ m_weakPhaseCoupling = val; }

        bool getViscosity(){ return m_viscosity; }
        void setViscosity(Real val){ m_viscosity = val; }

		void computeDensityGradient();
    };
}

#endif
