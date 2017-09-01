#ifndef __TimeStep_h__
#define __TimeStep_h__

#include "Common.h"
#include "FluidModel.h"
#include "NonPressureForceBase.h"
#include "Vorticity/MicropolarModel_Bender2017.h"
#include "Viscosity/ViscosityBase.h"
#include "SurfaceTension/SurfaceTensionBase.h"
#include "Vorticity/VorticityBase.h"
#include "Drag/DragBase.h"

namespace SPH
{
	enum class SurfaceTensionMethods { None = 0, Becker2007, Akinci2013, He2014 }; 
	enum class ViscosityMethods { None = 0, Standard, XSPH, Bender2017, Peer2015, Peer2016 };
	enum class VorticityMethods { None = 0, Micropolar, VorticityConfinement };
	enum class DragMethods { None = 0, Macklin2014, Gissler2017 };

	/** \brief Base class for the simulation methods. 
	*/
	class TimeStep
	{
	protected:
		FluidModel *m_model;
		unsigned int m_iterations;	
		unsigned int m_iterationsV;
		unsigned int m_cflMethod;
		Real m_cflFactor;
		Real m_cflMaxTimeStepSize;
		Real m_maxError;
		unsigned int m_maxIterations;
		Real m_maxErrorV;
		unsigned int m_maxIterationsV;
		SurfaceTensionMethods m_surfaceTensionMethod;
		SurfaceTensionBase *m_surfaceTension;
		ViscosityMethods m_viscosityMethod;
		ViscosityBase *m_viscosity;
		VorticityMethods m_fluidModel;
		VorticityBase *m_vorticity;
		DragMethods m_dragMethod;
		DragBase *m_drag;

		/** Clear accelerations and add gravitation.
		*/
		void clearAccelerations();

		/** Determine densities of all fluid particles.
		*/
		void computeDensities();

		/** Update time step size depending on the chosen method.
		*/
		void updateTimeStepSize();

		/** Update time step size by CFL condition.
		*/
		void updateTimeStepSizeCFL(const Real minTimeStepSize);

		/** Perform the neighborhood search for all fluid particles.
		*/
		virtual void performNeighborhoodSearch();
		void performNeighborhoodSearchSort();

		void computeSurfaceTension();
		void computeViscosity();
		void computeVorticity();
		void computeDragForce();
		void computeNonPressureForces();

	public:
		TimeStep(FluidModel *model);
		virtual ~TimeStep(void);

		virtual void step() = 0;
		virtual void reset();

		FluidModel *getModel() { return m_model; }
		unsigned int getIterationCount() const { return m_iterations; }
		unsigned int getIterationCountV() const { return m_iterationsV; }
		
		unsigned int getCflMethod() const { return m_cflMethod; }
		void setCflMethod(unsigned int val) { m_cflMethod = val; }
		Real getCflFactor() const { return m_cflFactor; }
		void setCflFactor(Real val) { m_cflFactor = val; }
		Real getCflMaxTimeStepSize() const { return m_cflMaxTimeStepSize; }
		void setCflMaxTimeStepSize(Real val) { m_cflMaxTimeStepSize = val; }
		Real getMaxError() const { return m_maxError; }
		void setMaxError(Real val) { m_maxError = val; }
		unsigned int getMaxIterations() const { return m_maxIterations; }
		void setMaxIterations(unsigned int val) { m_maxIterations = val; }
		Real getMaxErrorV() const { return m_maxErrorV; }
		void setMaxErrorV(Real val) { m_maxErrorV = val; }
		unsigned int getMaxIterationsV() const { return m_maxIterationsV; }
		void setMaxIterationsV(unsigned int val) { m_maxIterationsV = val; }
		SurfaceTensionMethods getSurfaceTensionMethod() const { return m_surfaceTensionMethod; }
		void setSurfaceTensionMethod(SurfaceTensionMethods val);
		ViscosityMethods getViscosityMethod() const { return m_viscosityMethod; }
		void setViscosityMethod(ViscosityMethods val);
		VorticityMethods getVorticityMethod() const { return m_fluidModel; }
		void setVorticityMethod(SPH::VorticityMethods val);
		DragMethods getDragMethod() const { return m_dragMethod; }
		void setDragMethod(SPH::DragMethods val);
		void emitParticles();

		virtual void emittedParticles(const unsigned int startIndex);

		SurfaceTensionBase *getSurfaceTensionBase() { return m_surfaceTension; }
		ViscosityBase *getViscosityBase() { return m_viscosity; }
		VorticityBase *getVorticityBase() { return m_vorticity; }
		DragBase *getDragBase() { return m_drag; }
	};
}

#endif
