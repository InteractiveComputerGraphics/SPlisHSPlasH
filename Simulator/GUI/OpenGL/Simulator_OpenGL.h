#ifndef __Simulator_OpenGL_h__
#define __Simulator_OpenGL_h__

#include "SPlisHSPlasH/Common.h"
#include "GUI/OpenGL/Shader.h"
#include <vector>

namespace SPH 
{	
	class FluidModel;
	class BoundaryModel;
	class BoundaryModel_Akinci2012;

	class Simulator_OpenGL
	{
	protected:
		static Shader m_shader_vector;
		static Shader m_shader_scalar;
		static Shader m_shader_scalar_map;
		static Shader m_meshShader;
		static GLuint m_textureMap;
	
	public:
		Simulator_OpenGL();
		virtual ~Simulator_OpenGL();
		
		static void initShaders(const std::string &shaderPath);
		static Shader& getShaderVector() { return m_shader_vector; }
		static Shader& getShaderScalar() { return m_shader_scalar; }
		static Shader& getMeshShader() { return m_meshShader; }
		static void meshShaderBegin(const float *col);
		static void meshShaderEnd();
		static void pointShaderBegin(Shader *shader, const Real particleRadius, const float *col, const Real minVal, const Real maxVal, const bool useTexture = false, float const* color_map = nullptr);
		static void pointShaderEnd(Shader *shader, const bool useTexture = false);
		static void renderFluid(FluidModel *model, float *fluidColor, 
			const unsigned int colorMapType, const bool useScalarField, const std::vector<float>& scalarField,
			const Real renderMinValue, const Real renderMaxValue);
		static void renderSelectedParticles(FluidModel *model, const std::vector<std::vector<unsigned int>>& selectedParticles,
			const unsigned int colorMapType, 
			const Real renderMinValue, const Real renderMaxValue);
		static void renderBoundary(BoundaryModel *model, const float *col);
		static void renderBoundaryParticles(BoundaryModel_Akinci2012 *model, const float *col);
	};
}

#endif