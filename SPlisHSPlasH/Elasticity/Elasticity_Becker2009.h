#ifndef __Elasticity_Becker2009_h__
#define __Elasticity_Becker2009_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ElasticityBase.h"

namespace SPH
{
	/** \brief This class implements the corotated SPH method for deformable solids introduced
	* by Becker et al. \cite Becker:2009.
	*/
	class Elasticity_Becker2009 : public ElasticityBase
	{
	protected:
		// initial particle indices, used to access their original positions
		std::vector<unsigned int> m_current_to_initial_index;
		std::vector<unsigned int> m_initial_to_current_index;
		// initial particle neighborhood
		std::vector<std::vector<unsigned int>> m_initialNeighbors;
		// volumes in rest configuration
		std::vector<Real> m_restVolumes;
		std::vector<Matrix3r> m_rotations;
		std::vector<Vector6r> m_stress;
		std::vector<Matrix3r> m_F;
		Real m_alpha;

		void initValues();
		void computeRotations();
		void computeStress();
		void computeForces();

		virtual void initParameters();

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
		static int ALPHA;

		Elasticity_Becker2009(FluidModel *model);
		virtual ~Elasticity_Becker2009(void);

		virtual void step();
		virtual void reset();
		virtual void performNeighborhoodSearchSort();

		virtual void saveState(BinaryFileWriter &binWriter);
		virtual void loadState(BinaryFileReader &binReader);
	};
}

#endif
