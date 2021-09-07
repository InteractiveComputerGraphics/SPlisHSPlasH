#pragma once

#include "SPHSamplingBase.h"
#include "SPlisHSPlasH/NeighborhoodSearch.h"


namespace SPH
{
	/** \brief This class implements the sampling method introduced
	* by Jiang et al. [JZW15].
	*
	* References:
	* - [JZW15] M. Jiang, Y. Zhou, R. Wang, R. Southern, J. J. Zhang.
	* Blue noise sampling using an SPH-based method. 
	* ACM Transactions on Graphics, 2015
	*/
	class SPHVolumeSampling_Jiang2015 : public SPHSamplingBase
	{
	protected:
		Real m_stiffness;
		Real m_dt;
		std::vector<Vector3r> m_v;
		std::vector<Vector3r> m_a;
		std::vector<Vector3r> m_n;
		std::vector<Real> m_p;

		void computePressure();
		void computeCohesion();
		void computeNormals();

		virtual void initSPHOptimization();
		virtual void step(Real &avg_density_error);	
		virtual void generateSamples();
		

	public:
		SPHVolumeSampling_Jiang2015();
		~SPHVolumeSampling_Jiang2015();

		Real getTimeStepSize() const { return m_dt; }
		void setTimeStepSize(const Real dt) { m_dt; }
		Real getStiffness() const { return m_stiffness; }
		void setStiffness(const Real stiffness) { m_stiffness = stiffness; }
	};
}
