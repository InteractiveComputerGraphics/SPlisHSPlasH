#ifndef __TimeStep_h__
#define __TimeStep_h__

#include "Common.h"
#include "ParameterObject.h"
#include "FluidModel.h"
#include "BoundaryModel.h"
#include "Discregrid/discrete_grid.hpp"

namespace SPH
{
	/** \brief Base class for the simulation methods. 
	*/
	class TimeStep : public GenParam::ParameterObject
	{
	public: 
		static int SOLVER_ITERATIONS;
		static int MIN_ITERATIONS;
		static int MAX_ITERATIONS;
		static int MAX_ERROR;

	protected:
		unsigned int m_iterations;	
		Real m_maxError;
		unsigned int m_minIterations;
		unsigned int m_maxIterations;	

		/** Clear accelerations and add gravitation.
		*/
		void clearAccelerations(const unsigned int fluidModelIndex);

		/** Determine densities of all fluid particles.
		*/
		void computeDensities(const unsigned int fluidModelIndex);

		virtual void initParameters();

		void approximateNormal(Discregrid::DiscreteGrid* map, const Eigen::Vector3d &x, Vector3r &n, const unsigned int dim);
		void computeVolumeAndBoundaryX(const unsigned int fluidModelIndex, const unsigned int i, const Vector3r &xi);
		void computeVolumeAndBoundaryX();
		void computeDensityAndGradient(const unsigned int fluidModelIndex, const unsigned int i, const Vector3r &xi);
		void computeDensityAndGradient();

	public:
		TimeStep();
		virtual ~TimeStep(void);

		virtual void step() = 0;
		virtual void reset();

		virtual void init();
		virtual void resize() = 0;

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex) {};

		virtual void saveState(BinaryFileWriter &binWriter) {};
		virtual void loadState(BinaryFileReader &binReader) {};
	};
}

#endif
