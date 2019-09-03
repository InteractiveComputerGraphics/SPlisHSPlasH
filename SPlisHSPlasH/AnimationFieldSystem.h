#ifndef __AnimationFieldSystem_h__
#define __AnimationFieldSystem_h__

#include "Common.h"
#include <vector>
#include "AnimationField.h"

namespace SPH 
{	
	class TimeStep;
	class FluidModel;

	class AnimationFieldSystem
	{
		public:
			AnimationFieldSystem();
			virtual ~AnimationFieldSystem();

		protected:
			std::vector<AnimationField*> m_fields;

			//void resetState();

		public:
			void addAnimationField(const std::string &particleFieldName, const Vector3r &pos, const Matrix3r & rotation, const Vector3r &scale,
				const std::string expression[3], const unsigned int type);
			unsigned int numAnimationFields() const { return static_cast<unsigned int>(m_fields.size()); }
			std::vector<AnimationField*> &getAnimationFields() { return m_fields; }

			void step();
			void reset();
	};
}

#endif