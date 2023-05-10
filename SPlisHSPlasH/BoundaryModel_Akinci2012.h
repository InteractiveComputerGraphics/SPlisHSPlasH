#ifndef __BoundaryModel_Akinci2012_h__
#define __BoundaryModel_Akinci2012_h__

#include "Common.h"
#include <vector>

#include "BoundaryModel.h"
#include "SPHKernels.h"
#include "FluidModel.h"


namespace SPH 
{	
	class TimeStep;

	/** \brief The boundary model stores the information required for boundary handling
	* using the approach of Akinci et al. 2012 [AIA+12].
	*
	* References:
	* - [AIA+12] Nadir Akinci, Markus Ihmsen, Gizem Akinci, Barbara Solenthaler, and Matthias Teschner. Versatile rigid-fluid coupling for incompressible SPH. ACM Trans. Graph., 31(4):62:1-62:8, July 2012. URL: http://doi.acm.org/10.1145/2185520.2185558
	*/
	class BoundaryModel_Akinci2012 : public BoundaryModel
	{
		public:
			BoundaryModel_Akinci2012();
			virtual ~BoundaryModel_Akinci2012();

		protected:
			bool m_sorted;
			unsigned int m_pointSetIndex;

			// values required for Akinci 2012 boundary handling
			std::vector<Vector3r> m_x0;
			std::vector<Vector3r> m_x;
			std::vector<Vector3r> m_v;
			std::vector<Real> m_V;
			std::vector<FieldDescription> m_fields;

		public:
			void addField(const FieldDescription& field);
			const std::vector<FieldDescription>& getFields() {return m_fields;}
			const FieldDescription& getField(const unsigned int i) {return m_fields[i];}
			const FieldDescription& getField(const std::string& name);
			const unsigned int numberOfFields() {return static_cast<unsigned int>(m_fields.size());}
			void removeFieldByName(const std::string& fieldName);

			unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }
			unsigned int getPointSetIndex() const { return m_pointSetIndex; }
			bool isSorted() const { return m_sorted; }

			void computeBoundaryVolume();
			void resize(const unsigned int numBoundaryParticles);

			virtual void reset();

			virtual void performNeighborhoodSearchSort();

			virtual void saveState(BinaryFileWriter &binWriter);
			virtual void loadState(BinaryFileReader &binReader);

			void initModel(RigidBodyObject *rbo, const unsigned int numBoundaryParticles, Vector3r *boundaryParticles);

		
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

			FORCE_INLINE const Real& getVolume(const unsigned int i) const
			{
				return m_V[i];
			}

			FORCE_INLINE Real& getVolume(const unsigned int i)
			{
				return m_V[i];
			}

			FORCE_INLINE void setVolume(const unsigned int i, const Real &val)
			{
				m_V[i] = val;
			}
	};
}

#endif