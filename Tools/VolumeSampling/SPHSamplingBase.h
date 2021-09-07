#pragma once

#include "SamplingBase.h"
#include "SPlisHSPlasH/NeighborhoodSearch.h"


namespace SPH
{
	class SPHSamplingBase : public SamplingBase
	{
	protected:
		Real m_viscosity;
		Real m_cohesion;
		Real m_adhesion;
		Real m_cflFactor;
		unsigned int m_steps;
		NeighborhoodSearch* m_neighborhoodSearch;

		std::vector<Real> m_densities;
		std::vector<Real> m_factors;
		std::vector<Vector3r> m_corrs;
		std::vector<Vector3r> m_deltaCorrs;

		const Real m_density0 = 1000.0;
		Real m_volume;
		Real m_mass;
		unsigned int m_counter;
		Real m_W_zero;
		Real(*m_kernelFct)(const Vector3r&);
		Vector3r(*m_gradKernelFct)(const Vector3r& r);

		virtual void initSPHOptimization() = 0;
		virtual void step(Real &avg_density_error) = 0;
		void writeParticleDataVTK(const std::string& fileName);

		void computeDFSPHFactor();
		void computeDensities(std::vector<Real>& densities, const Real mass);
		void pressureSolve();
		void pressureSolveIteration(Real& avg_density_err);
		
		virtual void generateSamples() = 0;


		FORCE_INLINE unsigned int numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int index)
		{
			return static_cast<unsigned int>(m_neighborhoodSearch->point_set(0).n_neighbors(pointSetIndex, index));
		}

		FORCE_INLINE unsigned int getNeighbor(const unsigned int pointSetIndex, const unsigned int index, const unsigned int k)
		{
			return m_neighborhoodSearch->point_set(0).neighbor(pointSetIndex, index, k);
		}

	public:
		SPHSamplingBase();
		~SPHSamplingBase();

		Real getViscosity() const { return m_viscosity;	}
		void setViscosity(const Real viscosity)	{ m_viscosity = viscosity; }
		Real getCohesion() const { return m_cohesion; }
		void setCohesion(const Real cohesion) { m_cohesion = cohesion; }
		Real getAdhesion() const { return m_adhesion; }
		void setAdhesion(const Real adhesion) { m_adhesion = adhesion; }
		Real getCFLFactor() const { return m_cflFactor; }
		void setCFLFactor(const Real cfl_factor) { m_cflFactor = cfl_factor; }
		unsigned getSteps() const { return m_steps; }
		void setSteps(const unsigned steps) { m_steps = steps; }
	};
}
