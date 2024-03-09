#ifndef __PBD_Simulator_GUI_imgui_h__
#define __PBD_Simulator_GUI_imgui_h__

#include "SPlisHSPlasH/Common.h"
#include "Simulator/GUI/imgui/Simulator_GUI_imgui.h"
#include "PBDWrapper.h"
#include "GUI/OpenGL/Shader.h"
#include "GUI/OpenGL/MiniGL.h"

namespace SPH 
{	
	class PBD_Simulator_GUI_imgui : public Simulator_GUI_imgui
	{
		protected: 
			PBDWrapper *m_pbdWrapper;
			bool m_drawAABB;
			int m_drawBVHDepth;
			bool m_drawSDF;
			SPH::Shader *m_shader;
			float m_jointColor[4];

			SPH::Shader *createShader(const std::string &vertexShader, const std::string &geometryShader, const std::string &fragmentShader);
			void initShader();

			void shaderBegin(const float *col);
			void shaderEnd();

			void renderAABB(PBD::AABB &aabb);
			void renderBVH();
			void renderSDF();
			void renderSDF(PBD::CollisionDetection::CollisionObject* co);
			void renderTriangleModels();
			void renderTetModels();
			void renderConstraints();
			void renderBallJoint(PBD::BallJoint &bj);
			void renderRigidBodyParticleBallJoint(PBD::RigidBodyParticleBallJoint &bj);
			void renderBallOnLineJoint(PBD::BallOnLineJoint &bj);
			void renderHingeJoint(PBD::HingeJoint &hj);
			void renderUniversalJoint(PBD::UniversalJoint &uj);
			void renderSliderJoint(PBD::SliderJoint &joint);
			void renderTargetPositionMotorSliderJoint(PBD::TargetPositionMotorSliderJoint &joint);
			void renderTargetVelocityMotorSliderJoint(PBD::TargetVelocityMotorSliderJoint &joint);
			void renderTargetAngleMotorHingeJoint(PBD::TargetAngleMotorHingeJoint &hj);
			void renderTargetVelocityMotorHingeJoint(PBD::TargetVelocityMotorHingeJoint &hj);
			void renderRigidBodyContact(PBD::RigidBodyContactConstraint &cc);
			void renderParticleRigidBodyContact(PBD::ParticleRigidBodyContactConstraint &cc);
			void renderSpring(PBD::RigidBodySpring &s);
			void renderDistanceJoint(PBD::DistanceJoint &j);
			void renderDamperJoint(PBD::DamperJoint &joint);

		public:
			PBD_Simulator_GUI_imgui(SimulatorBase *base, PBDWrapper *pbdWrapper);
			virtual ~PBD_Simulator_GUI_imgui();

		public:
			virtual void init(const char *name);
			virtual void render();
			virtual void initSimulationParameterGUI();

		public:
			template<class PositionData>
			static void drawMesh(const PositionData &pd, const Utilities::IndexedFaceMesh &mesh, const unsigned int offset, const float * const color);

	};


	template<class PositionData>
	void PBD_Simulator_GUI_imgui::drawMesh(const PositionData &pd, const Utilities::IndexedFaceMesh &mesh, const unsigned int offset, const float * const color)
	{
		// draw mesh 
		const unsigned int *faces = mesh.getFaces().data();
		const unsigned int nFaces = mesh.numFaces();
		const Vector3r *vertexNormals = mesh.getVertexNormals().data();

		SPH::MiniGL::supplyVertices(0, mesh.getVertexNormals().size(), &pd.getPosition(offset)[0]);
		SPH::MiniGL::supplyNormals(2, mesh.getVertexNormals().size(), &vertexNormals[0][0]);
		SPH::MiniGL::supplyFaces(3 * nFaces, faces);

		glDrawElements(GL_TRIANGLES, (GLsizei)3 * nFaces, GL_UNSIGNED_INT, (void*)0);

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(2);
	}

}

#endif