#ifndef __AnimationField_h__
#define __AnimationField_h__

#include "Common.h"
#include <vector>
#include "FluidModel.h"


namespace SPH 
{	
	class AnimationField 
	{
		public:
			AnimationField(
				const std::string &particleFieldName,
				const Vector3r &pos, const Matrix3r & rotation, const Vector3r &scale,
				const std::string expression[3], const unsigned int type = 0);
			virtual ~AnimationField();

		protected:
			std::string m_particleFieldName;
			Vector3r m_x;
			Matrix3r m_rotation;
			Vector3r m_scale;
			std::string m_expression[3];
			unsigned int m_type;
			Real m_startTime;
			Real m_endTime;

			FORCE_INLINE bool inBox(const Vector3r &x, const Vector3r &xBox, const Matrix3r &rotBox, const Vector3r &scaleBox)
			{
				const Vector3r xlocal = rotBox.transpose() * (x - xBox);
				// for a box shape, m_scale stores the half-size of the box
				// inside box if closer than half-size on all axes
				return (xlocal.array().abs() < scaleBox.array()).all();
			}

			FORCE_INLINE bool inCylinder(const Vector3r &x, const Vector3r &xCyl, const Matrix3r &rotCyl, const Real h, const Real r2)
			{
				const Vector3r xlocal = rotCyl.transpose() * (x - xCyl);
				// inside cylinder if distance to x-axis is less than r
				// and projection on x-axis is between 0 and h
				const Real proj = xlocal.x();
				const Real d2 = Vector2r(xlocal.y(), xlocal.z()).squaredNorm();
				const Real hHalf = 0.5*h;
				return (proj > -hHalf) && (proj < hHalf) && (d2 < r2);
			}

			FORCE_INLINE bool inSphere(const Vector3r &x, const Vector3r &pos, const Matrix3r &rot, const Real radius)
			{
				const Vector3r xlocal = rot.transpose() * (x - pos);
				return (xlocal.norm() < radius);
			}

			FORCE_INLINE bool inShape(const int type, const Vector3r &x, const Vector3r &pos, const Matrix3r &rot, const Vector3r &scale)
			{
				if (type == 1)
					return inSphere(x, pos, rot, scale[0]);
				else if (type == 2)
				{
					const Real h = scale[0];
					const Real r = scale[1];
					return inCylinder(x, pos, rot, h, r*r);
				}
				else 
					return inBox(x, pos, rot, 0.5*scale);
			}

		public:
			void setStartTime(Real val) { m_startTime = val; }
			void setEndTime(Real val) { m_endTime = val; }

			void step();
			virtual void reset();

	};
}

#endif