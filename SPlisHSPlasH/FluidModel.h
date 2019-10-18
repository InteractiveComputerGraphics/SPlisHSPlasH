#ifndef __FluidModel_h__
#define __FluidModel_h__

#include "Common.h"
#include <vector>

#include "RigidBodyObject.h"
#include "SPHKernels.h"
#include "ParameterObject.h"
#include "Utilities/BinaryFileReaderWriter.h"

namespace SPH 
{	
	class TimeStep;
	class ViscosityBase;
	class SurfaceTensionBase;
	class VorticityBase;
	class DragBase;
	class ElasticityBase;
	class EmitterSystem;

	enum FieldType { Scalar = 0, Vector3, Vector6, Matrix3, Matrix6, UInt };
	struct FieldDescription
	{
		std::string name;
		FieldType type;
		// getFct(particleIndex)
		std::function<void*(const unsigned int)> getFct;
		bool storeData;

		FieldDescription(const std::string &n, const FieldType &t, 
			const std::function<void*(const unsigned int)> &fct, const bool s = false) :
			name(n), type(t), getFct(fct), storeData(s) { }
	};

	enum class SurfaceTensionMethods { None = 0, Becker2007, Akinci2013, He2014, NumSurfaceTensionMethods };
	enum class ViscosityMethods { None = 0, Standard, XSPH, Bender2017, Peer2015, Peer2016, Takahashi2015, Weiler2018, NumViscosityMethods };
	enum class VorticityMethods { None = 0, Micropolar, VorticityConfinement, NumVorticityMethods };
	enum class DragMethods { None = 0, Macklin2014, Gissler2017, NumDragMethods };
	enum class ElasticityMethods { None = 0, Becker2009, Peer2018, NumElasticityMethods };

	enum class ParticleState { Active = 0, AnimatedByEmitter };

	/** \brief The fluid model stores the particle and simulation information 
	*/
	class FluidModel : public GenParam::ParameterObject
	{
		public:
			static int NUM_PARTICLES;
			static int NUM_REUSED_PARTICLES;
			static int DENSITY0;

			static int DRAG_METHOD;
			static int SURFACE_TENSION_METHOD;
			static int VISCOSITY_METHOD;
			static int VORTICITY_METHOD;
			static int ELASTICITY_METHOD;

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

			static int ENUM_ELASTICITY_NONE;
			static int ENUM_ELASTICITY_BECKER2009;
			static int ENUM_ELASTICITY_PEER2018;
			
			FluidModel();
			virtual ~FluidModel();

			void init();

			std::string getId() const { return m_id; }

		protected:
			std::string m_id;
			EmitterSystem *m_emitterSystem;

			// Mass
			// If the mass is zero, the particle is static
			std::vector<Real> m_masses;
			std::vector<Vector3r> m_a;
			std::vector<Vector3r> m_v0;
			std::vector<Vector3r> m_x0;
			std::vector<Vector3r> m_x;
			std::vector<Vector3r> m_v;
			std::vector<Real> m_density;
			std::vector<unsigned int> m_particleId;
			std::vector<ParticleState> m_particleState;
			Real m_V;

			SurfaceTensionMethods m_surfaceTensionMethod;
			SurfaceTensionBase *m_surfaceTension;
			ViscosityMethods m_viscosityMethod;
			ViscosityBase *m_viscosity;
			VorticityMethods m_vorticityMethod;
			VorticityBase *m_vorticity;
			DragMethods m_dragMethod;
			DragBase *m_drag;
			ElasticityMethods m_elasticityMethod;
			ElasticityBase *m_elasticity;
			std::vector<FieldDescription> m_fields;

			std::function<void()> m_dragMethodChanged;
			std::function<void()> m_surfaceTensionMethodChanged;
			std::function<void()> m_viscosityMethodChanged;
			std::function<void()> m_vorticityMethodChanged;
			std::function<void()> m_elasticityMethodChanged;

			Real m_density0;
			unsigned int m_pointSetIndex;

			unsigned int m_numActiveParticles;
			unsigned int m_numActiveParticles0;

			virtual void initParameters();

			void initMasses();

			/** Resize the arrays containing the particle data.
			*/
			virtual void resizeFluidParticles(const unsigned int newSize);
			
			/** Release the arrays containing the particle data.
			*/
			virtual void releaseFluidParticles();


		public:
			FORCE_INLINE Real getDensity0() const { return m_density0; }
			void setDensity0(const Real v);

			unsigned int getPointSetIndex() const { return m_pointSetIndex; }

			void addField(const FieldDescription &field);
			const std::vector<FieldDescription> &getFields() { return m_fields; }
			const FieldDescription &getField(const unsigned int i) { return m_fields[i]; }
			const FieldDescription &getField(const std::string &name);
			const unsigned int numberOfFields() { return static_cast<unsigned int>(m_fields.size()); }
			void removeFieldByName(const std::string &fieldName);

			void setNumActiveParticles(const unsigned int num);
			unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }

			EmitterSystem* getEmitterSystem() { return m_emitterSystem; }

			virtual void reset();

			void performNeighborhoodSearchSort();

			void initModel(const std::string &id, const unsigned int nFluidParticles, Vector3r* fluidParticles, Vector3r* fluidVelocities, const unsigned int nMaxEmitterParticles);
			
			const unsigned int numParticles() const { return static_cast<unsigned int>(m_masses.size()); }
			unsigned int numActiveParticles() const;

			unsigned int getNumActiveParticles0() const { return m_numActiveParticles0; }
			void setNumActiveParticles0(unsigned int val) { m_numActiveParticles0 = val; }

			void emittedParticles(const unsigned int startIndex);

			int getSurfaceTensionMethod() const { return static_cast<int>(m_surfaceTensionMethod); }
			void setSurfaceTensionMethod(const int val);
			int getViscosityMethod() const { return static_cast<int>(m_viscosityMethod); }
			void setViscosityMethod(const int val);
			int getVorticityMethod() const { return static_cast<int>(m_vorticityMethod); }
			void setVorticityMethod(const int val);
			int getDragMethod() const { return static_cast<int>(m_dragMethod); }
			void setDragMethod(const int val);
			int getElasticityMethod() const { return static_cast<int>(m_elasticityMethod); }
			void setElasticityMethod(const int val);

			SurfaceTensionBase *getSurfaceTensionBase() { return m_surfaceTension; }
			ViscosityBase *getViscosityBase() { return m_viscosity; }
			VorticityBase *getVorticityBase() { return m_vorticity; }
			DragBase *getDragBase() { return m_drag; }
			ElasticityBase *getElasticityBase() { return m_elasticity; }

			void setDragMethodChangedCallback(std::function<void()> const& callBackFct);
			void setSurfaceMethodChangedCallback(std::function<void()> const& callBackFct);
			void setViscosityMethodChangedCallback(std::function<void()> const& callBackFct);
			void setVorticityMethodChangedCallback(std::function<void()> const& callBackFct);
			void setElasticityMethodChangedCallback(std::function<void()> const& callBackFct);

			void computeSurfaceTension();
			void computeViscosity();
			void computeVorticity();
			void computeDragForce();
			void computeElasticity();

			void saveState(BinaryFileWriter &binWriter);
			void loadState(BinaryFileReader &binReader);


			FORCE_INLINE Vector3r &getPosition0(const unsigned int i)
			{
				return m_x0[i];
			}

			FORCE_INLINE const Vector3r &getPosition0(const unsigned int i) const
			{
				return m_x0[i];
			}

			FORCE_INLINE void setPosition0(const unsigned int i, const Vector3r &pos)
			{
				m_x0[i] = pos;
			}

			FORCE_INLINE Vector3r &getPosition(const unsigned int i)
			{
				return m_x[i];
			}

			FORCE_INLINE const Vector3r &getPosition(const unsigned int i) const
			{
				return m_x[i];
			}

			FORCE_INLINE void setPosition(const unsigned int i, const Vector3r &pos)
			{
				m_x[i] = pos;
			}

			FORCE_INLINE Vector3r &getVelocity(const unsigned int i)
			{
				return m_v[i];
			}

			FORCE_INLINE const Vector3r &getVelocity(const unsigned int i) const
			{
				return m_v[i];
			}

			FORCE_INLINE void setVelocity(const unsigned int i, const Vector3r &vel)
			{
				m_v[i] = vel;
			}

			FORCE_INLINE Vector3r &getVelocity0(const unsigned int i)
			{
				return m_v0[i];
			}

			FORCE_INLINE const Vector3r &getVelocity0(const unsigned int i) const
			{
				return m_v0[i];
			}

			FORCE_INLINE void setVelocity0(const unsigned int i, const Vector3r &vel)
			{
				m_v0[i] = vel;
			}

			FORCE_INLINE Vector3r &getAcceleration(const unsigned int i)
			{
				return m_a[i];
			}

			FORCE_INLINE const Vector3r &getAcceleration(const unsigned int i) const
			{
				return m_a[i];
			}

			FORCE_INLINE void setAcceleration(const unsigned int i, const Vector3r &accel)
			{
				m_a[i] = accel;
			}

			FORCE_INLINE const Real getMass(const unsigned int i) const
			{
				return m_masses[i];
			}

			FORCE_INLINE Real& getMass(const unsigned int i)
			{
				return m_masses[i];
			}

			FORCE_INLINE void setMass(const unsigned int i, const Real mass)
			{
				m_masses[i] = mass;
			}

			FORCE_INLINE const Real& getDensity(const unsigned int i) const
			{
				return m_density[i];
			}

			FORCE_INLINE Real& getDensity(const unsigned int i)
			{
				return m_density[i];
			}

			FORCE_INLINE void setDensity(const unsigned int i, const Real &val)
			{
				m_density[i] = val;
			}

			FORCE_INLINE unsigned int& getParticleId(const unsigned int i)
			{
				return m_particleId[i];
			}

			FORCE_INLINE const unsigned int& getParticleId(const unsigned int i) const
			{
				return m_particleId[i];
			}

			FORCE_INLINE const ParticleState& getParticleState(const unsigned int i) const
			{
				return m_particleState[i];
			}

			FORCE_INLINE ParticleState& getParticleState(const unsigned int i)
			{
				return m_particleState[i];
			}

			FORCE_INLINE void setParticleState(const unsigned int i, const ParticleState &val)
			{
				m_particleState[i] = val;
			}

			FORCE_INLINE const Real getVolume(const unsigned int i) const
			{
				return m_V;
			}

			FORCE_INLINE Real& getVolume(const unsigned int i)
			{
				return m_V;
			}
	};
}

#endif