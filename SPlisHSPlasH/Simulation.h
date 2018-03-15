#ifndef __Simulation_h__
#define __Simulation_h__

#include "Common.h"
#include "Viscosity/ViscosityBase.h"
#include "SurfaceTension/SurfaceTensionBase.h"
#include "Vorticity/VorticityBase.h"
#include "Drag/DragBase.h"
#include "FluidModel.h"
#include "NonPressureForceBase.h"
#include "ParameterObject.h"


namespace SPH
{
	enum class SimulationMethods { WCSPH = 0, PCISPH, PBF, IISPH, DFSPH, PF, NumSimulationMethods };
	enum class SurfaceTensionMethods { None = 0, Becker2007, Akinci2013, He2014, NumSurfaceTensionMethods };
	enum class ViscosityMethods { None = 0, Standard, XSPH, Bender2017, Peer2015, Peer2016, Takahashi2015, Weiler2018, NumViscosityMethods };
	enum class VorticityMethods { None = 0, Micropolar, VorticityConfinement, NumVorticityMethods };
	enum class DragMethods { None = 0, Macklin2014, Gissler2017, NumDragMethods };

	/** \brief Class to manage the current simulation time and the time step size. 
	* This class is a singleton.
	*/
	class Simulation : public GenParam::ParameterObject
	{
	public:
		static int GRAVITATION;
		static int CFL_METHOD;
		static int CFL_FACTOR;
		static int CFL_MAX_TIMESTEPSIZE;
		static int DRAG_METHOD;
		static int SURFACE_TENSION_METHOD;
		static int VISCOSITY_METHOD;
		static int VORTICITY_METHOD;
		static int SIMULATION_METHOD;

		static int ENUM_CFL_NONE;
		static int ENUM_CFL_STANDARD;
		static int ENUM_CFL_ITER;

		static int ENUM_DRAG_NONE;
		static int ENUM_DRAG_MACKLIN2014;
		static int ENUM_DRAG_GISSLER2017;

		static int ENUM_SURFACETENSION_NONE;
		static int ENUM_SURFACETENSION_BECKER2007;
		static int ENUM_SURFACETENSION_AKINCI2013;
		static int ENUM_SURFACETENSION_HE2014;

		static int ENUM_VISCOSITY_NONE;
		static int ENUM_VISCOSITY_STANDARD;
		static int ENUM_VISCOSITY_XSPH;
		static int ENUM_VISCOSITY_BENDER2017;
		static int ENUM_VISCOSITY_PEER2015;
		static int ENUM_VISCOSITY_PEER2016;
		static int ENUM_VISCOSITY_TAKAHASHI2015;
		static int ENUM_VISCOSITY_WEILER2018;

		static int ENUM_VORTICITY_NONE;
		static int ENUM_VORTICITY_MICROPOLAR;
		static int ENUM_VORTICITY_VC;

		static int ENUM_SIMULATION_WCSPH;
		static int ENUM_SIMULATION_PCISPH;
		static int ENUM_SIMULATION_PBF;
		static int ENUM_SIMULATION_IISPH;
		static int ENUM_SIMULATION_DFSPH;
		static int ENUM_SIMULATION_PF;

	protected:
		FluidModel *m_model;
		int m_cflMethod;
		Real m_cflFactor;
		Real m_cflMaxTimeStepSize;
		SurfaceTensionMethods m_surfaceTensionMethod;
		SurfaceTensionBase *m_surfaceTension;
		ViscosityMethods m_viscosityMethod;
		ViscosityBase *m_viscosity;
		VorticityMethods m_vorticityMethod;
		VorticityBase *m_vorticity;
		DragMethods m_dragMethod;
		DragBase *m_drag;
		SimulationMethods m_simulationMethod;
		TimeStep *m_timeStep;
		Vector3r m_gravitation;
		std::function<void()> m_dragMethodChanged;
		std::function<void()> m_surfaceTensionMethodChanged;
		std::function<void()> m_viscosityMethodChanged;
		std::function<void()> m_vorticityMethodChanged;
		std::function<void()> m_simulationMethodChanged;

		virtual void initParameters();
		

	private:
		static Simulation *current;

	public:
		Simulation ();
		~Simulation ();

		void init();
		void reset();

		// Singleton
		static Simulation* getCurrent ();
		static void setCurrent (Simulation* tm);
		static bool hasCurrent();

		FluidModel *getModel() { return m_model; }

		int getSimulationMethod() const { return static_cast<int>(m_simulationMethod); }
		void setSimulationMethod(const int val);
		int getSurfaceTensionMethod() const { return static_cast<int>(m_surfaceTensionMethod); }
		void setSurfaceTensionMethod(const int val);
		int getViscosityMethod() const { return static_cast<int>(m_viscosityMethod); }
		void setViscosityMethod(const int val);
		int getVorticityMethod() const { return static_cast<int>(m_vorticityMethod); }
		void setVorticityMethod(const int val);
		int getDragMethod() const { return static_cast<int>(m_dragMethod); }
		void setDragMethod(const int val);

		void setDragMethodChangedCallback(std::function<void()> const& callBackFct);
		void setSurfaceMethodChangedCallback(std::function<void()> const& callBackFct);
		void setViscosityMethodChangedCallback(std::function<void()> const& callBackFct);
		void setVorticityMethodChangedCallback(std::function<void()> const& callBackFct);
		void setSimulationMethodChangedCallback(std::function<void()> const& callBackFct);

		SurfaceTensionBase *getSurfaceTensionBase() { return m_surfaceTension; }
		ViscosityBase *getViscosityBase() { return m_viscosity; }
		VorticityBase *getVorticityBase() { return m_vorticity; }
		DragBase *getDragBase() { return m_drag; }
		TimeStep *getTimeStep() { return m_timeStep; }

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

		void emitParticles();
		virtual void emittedParticles(const unsigned int startIndex);
	};
}

#endif
