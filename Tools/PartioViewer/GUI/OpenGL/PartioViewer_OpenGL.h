#pragma once

#include "GUI/OpenGL/Shader.h"
#include "SPlisHSPlasH/Common.h"
#include <vector>

namespace SPH
{
	struct Fluid;
	struct Boundary;

	class PartioViewer_OpenGL
	{
	protected:
		static Shader m_shader_vector;
		static Shader m_shader_scalar;
		static Shader m_shader_scalar_map;
		static Shader m_meshShader;
		static GLuint m_textureMap;

		static bool m_updateScalarField;
		static std::vector<float> m_scalarField;

	public:
		PartioViewer_OpenGL();
		~PartioViewer_OpenGL();

		static void flipImage(int width, int height, unsigned char *image);
		static void getImage(int width, int height, unsigned char *image);
		static void initShaders(const std::string &shaderPath);
		static void pointShaderBegin(Shader *shader, const Real particleRadius, const float *col, const Real minVal, const Real maxVal, const bool useTexture = false, float const* color_map = nullptr);
		static void pointShaderEnd(Shader *shader, const bool useTexture = false);
		static void renderGrid();
		static void renderAABB(const Eigen::AlignedBox3f &aabb, float *color);
		static void render(const Fluid &fluid, const Real particleRadius, 
			float *fluidColor, const unsigned int colorMapType, const unsigned int colorField, 
			const float renderMinValue, const float renderMaxValue, const bool usePlane);
		static void render(const Boundary &boundary, const bool renderWalls, float *boundaryColor);
		static void hsvToRgb(float h, float s, float v, float *rgb);
		static void updateScalarField();
	};
}
