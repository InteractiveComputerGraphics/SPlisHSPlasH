#include "PartioViewer_OpenGL.h"

#include "GUI/OpenGL/MiniGL.h"
#include "GUI/OpenGL/colormaps/colormap_jet.h"
#include "GUI/OpenGL/colormaps/colormap_plasma.h"
#include "GUI/OpenGL/colormaps/colormap_bwr.h"
#include "GUI/OpenGL/colormaps/colormap_coolwarm.h"
#include "GUI/OpenGL/colormaps/colormap_seismic.h"
#include "../../PartioViewer.h"

using namespace std;
using namespace SPH;

Shader PartioViewer_OpenGL::m_shader_scalar;
Shader PartioViewer_OpenGL::m_shader_scalar_map;
Shader PartioViewer_OpenGL::m_shader_vector;
Shader PartioViewer_OpenGL::m_meshShader;
GLuint PartioViewer_OpenGL::m_textureMap = 0;
std::vector<float> PartioViewer_OpenGL::m_scalarField;
bool PartioViewer_OpenGL::m_updateScalarField = true;

PartioViewer_OpenGL::PartioViewer_OpenGL()
{
}


PartioViewer_OpenGL::~PartioViewer_OpenGL()
{
}

void PartioViewer_OpenGL::flipImage(int width, int height, unsigned char *image)
{
	unsigned char tmp[3];

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height / 2; y++)
		{
			int top = (x + y * width) * 3;
			int bottom = (x + (height - y - 1) * width) * 3;

			memcpy(tmp, image + top, sizeof(tmp));
			memcpy(image + top, image + bottom, sizeof(tmp));
			memcpy(image + bottom, tmp, sizeof(tmp));
		}
	}
}

void PartioViewer_OpenGL::getImage(int width, int height, unsigned char *image)
{
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadBuffer(GL_FRONT_LEFT);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
	flipImage(width, height, image);
}


void PartioViewer_OpenGL::initShaders(const std::string &shaderPath)
{
	string vertFile = shaderPath + "/vs_points_vector.glsl";
	string fragFile = shaderPath + "/fs_points.glsl";
	m_shader_vector.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	m_shader_vector.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	m_shader_vector.createAndLinkProgram();
	m_shader_vector.begin();
	m_shader_vector.addUniform("modelview_matrix");
	m_shader_vector.addUniform("projection_matrix");
	m_shader_vector.addUniform("radius");
	m_shader_vector.addUniform("viewport_width");
	m_shader_vector.addUniform("color");
	m_shader_vector.addUniform("min_scalar");
	m_shader_vector.addUniform("max_scalar");
	m_shader_vector.end();

	string vertFileScalar = shaderPath + "/vs_points_scalar.glsl";
	m_shader_scalar.compileShaderFile(GL_VERTEX_SHADER, vertFileScalar);
	m_shader_scalar.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	m_shader_scalar.createAndLinkProgram();
	m_shader_scalar.begin();
	m_shader_scalar.addUniform("modelview_matrix");
	m_shader_scalar.addUniform("projection_matrix");
	m_shader_scalar.addUniform("radius");
	m_shader_scalar.addUniform("viewport_width");
	m_shader_scalar.addUniform("color");
	m_shader_scalar.addUniform("min_scalar");
	m_shader_scalar.addUniform("max_scalar");
	m_shader_scalar.end();

	string fragFileMap = shaderPath + "/fs_points_colormap.glsl";
	m_shader_scalar_map.compileShaderFile(GL_VERTEX_SHADER, vertFileScalar);
	m_shader_scalar_map.compileShaderFile(GL_FRAGMENT_SHADER, fragFileMap);
	m_shader_scalar_map.createAndLinkProgram();
	m_shader_scalar_map.begin();
	m_shader_scalar_map.addUniform("modelview_matrix");
	m_shader_scalar_map.addUniform("projection_matrix");
	m_shader_scalar_map.addUniform("radius");
	m_shader_scalar_map.addUniform("viewport_width");
	m_shader_scalar_map.addUniform("color");
	m_shader_scalar_map.addUniform("min_scalar");
	m_shader_scalar_map.addUniform("max_scalar");
	m_shader_scalar_map.end();

	vertFile = shaderPath + "/vs_smooth.glsl";
	fragFile = shaderPath + "/fs_smooth.glsl";
	m_meshShader.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	m_meshShader.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	m_meshShader.createAndLinkProgram();
	m_meshShader.begin();
	m_meshShader.addUniform("modelview_matrix");
	m_meshShader.addUniform("projection_matrix");
	m_meshShader.addUniform("surface_color");
	m_meshShader.addUniform("shininess");
	m_meshShader.addUniform("specular_factor");
	m_meshShader.end();

	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &m_textureMap);
}


void PartioViewer_OpenGL::pointShaderBegin(Shader *shader, const Real particleRadius, const float *col, const Real minVal, const Real maxVal, const bool useTexture, float const* color_map)
{
	shader->begin();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glUniform1f(shader->getUniform("viewport_width"), (float)viewport[2]);
	glUniform1f(shader->getUniform("radius"), (float)particleRadius);
	glUniform1f(shader->getUniform("min_scalar"), (GLfloat)minVal);
	glUniform1f(shader->getUniform("max_scalar"), (GLfloat)maxVal);
	glUniform3fv(shader->getUniform("color"), 1, col);

	if (useTexture)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_1D, m_textureMap);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 256u, 0, GL_RGB, GL_FLOAT, color_map);

		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glGenerateMipmap(GL_TEXTURE_1D);
	}


	GLfloat matrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
	glUniformMatrix4fv(shader->getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
	GLfloat pmatrix[16];
	glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
	glUniformMatrix4fv(shader->getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

	glEnable(GL_DEPTH_TEST);
	// Point sprites do not have to be explicitly enabled since OpenGL 3.2 where
	// they are enabled by default. Moreover GL_POINT_SPRITE is deprecate and only
	// supported before OpenGL 3.2 or with compatibility profile enabled.
	glEnable(GL_POINT_SPRITE);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glPointParameterf(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);
}

void PartioViewer_OpenGL::pointShaderEnd(Shader *shader, const bool useTexture)
{
	glBindTexture(GL_TEXTURE_1D, 0);
	shader->end();
}

void PartioViewer_OpenGL::renderAABB(const Eigen::AlignedBox3f &aabb, float *color)
{
	const Vector3r a = aabb.corner(Eigen::AlignedBox3f::BottomLeftFloor).cast<Real>();
	const Vector3r b = aabb.corner(Eigen::AlignedBox3f::BottomRightFloor).cast<Real>();
	const Vector3r c = aabb.corner(Eigen::AlignedBox3f::TopRightFloor).cast<Real>();
	const Vector3r d = aabb.corner(Eigen::AlignedBox3f::TopLeftFloor).cast<Real>();
	const Vector3r e = aabb.corner(Eigen::AlignedBox3f::BottomLeftCeil).cast<Real>();
	const Vector3r f = aabb.corner(Eigen::AlignedBox3f::BottomRightCeil).cast<Real>();
	const Vector3r g = aabb.corner(Eigen::AlignedBox3f::TopRightCeil).cast<Real>();
	const Vector3r h = aabb.corner(Eigen::AlignedBox3f::TopLeftCeil).cast<Real>();

	const float w = 1.0;
	MiniGL::drawVector(a, b, w, color);
	MiniGL::drawVector(b, c, w, color);
	MiniGL::drawVector(c, d, w, color);
	MiniGL::drawVector(d, a, w, color);

	MiniGL::drawVector(e, f, w, color);
	MiniGL::drawVector(f, g, w, color);
	MiniGL::drawVector(g, h, w, color);
	MiniGL::drawVector(h, e, w, color);

	MiniGL::drawVector(a, e, w, color);
	MiniGL::drawVector(b, f, w, color);
	MiniGL::drawVector(c, g, w, color);
	MiniGL::drawVector(d, h, w, color);
}

void PartioViewer_OpenGL::renderGrid()
{
	MiniGL::coordinateSystem();

	float gridColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
	MiniGL::drawGrid_xz(gridColor);
}

void PartioViewer_OpenGL::render(const Fluid &fluid, const Real particleRadius, 
		float *fluidColor, const unsigned int colorMapType, const unsigned int colorField,
		const float renderMinValue, const float renderMaxValue, const bool usePlane)
{
	// Draw simulation model
	const unsigned int nParticles = (unsigned int)fluid.partioData->numParticles();

	Partio::ParticleAttribute posAttr;
	fluid.partioData->attributeInfo(fluid.posIndex, posAttr);
	const float* partioX = fluid.partioData->data<float>(posAttr, 0);

	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		Shader *shader_s = &m_shader_scalar_map;
		float const *color_map = nullptr;
		if (colorMapType == 1)
			color_map = reinterpret_cast<float const*>(colormap_jet);
		else if (colorMapType == 2)
			color_map = reinterpret_cast<float const*>(colormap_plasma);
		else if (colorMapType == 3)
			color_map = reinterpret_cast<float const*>(colormap_coolwarm);
		else if (colorMapType == 4)
			color_map = reinterpret_cast<float const*>(colormap_bwr);
		else if (colorMapType == 5)
			color_map = reinterpret_cast<float const*>(colormap_seismic);

		if (colorMapType == 0)
			shader_s = &m_shader_scalar;

		if (fluid.partioData->numAttributes() == 0)
			pointShaderBegin(shader_s, particleRadius, fluidColor, renderMinValue, renderMaxValue, false);
		else
			pointShaderBegin(shader_s, particleRadius, fluidColor, renderMinValue, renderMaxValue, true, color_map);

		if (nParticles > 0)
		{
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, partioX);

			if (fluid.partioData->numAttributes() > 0)
			{
				glEnableVertexAttribArray(1);

				if (m_updateScalarField)
				{
					Partio::ParticleAttribute attr;
					fluid.partioData->attributeInfo(colorField, attr);

					m_scalarField.resize(nParticles);
					for (unsigned int i = 0u; i < nParticles; i++)
					{
						if (attr.type == Partio::VECTOR)
						{
							const Eigen::Map<const Eigen::Vector3f> vec(fluid.partioData->data<float>(attr, i));
							m_scalarField[i] = static_cast<float>(vec.norm());
						}
						else if (attr.type == Partio::FLOAT)
						{
							m_scalarField[i] = *fluid.partioData->data<float>(attr, i);
						}
						else if (attr.type == Partio::INT)
						{
							m_scalarField[i] = static_cast<float>(*fluid.partioData->data<int>(attr, i));
						}
					}
					m_updateScalarField = false;
					glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, &m_scalarField[0]);
				}
			}

			if (usePlane)
				glDrawElements(GL_POINTS, (GLsizei)fluid.visibleParticles.size(), GL_UNSIGNED_INT, fluid.visibleParticles.data());
			else
				glDrawArrays(GL_POINTS, 0, nParticles);

			glDisableVertexAttribArray(0);
			//			glDisableVertexAttribArray(1);
		}

		if (fluid.partioData->numAttributes() == 0)
			pointShaderEnd(shader_s, false);
		else
			pointShaderEnd(shader_s, true);
	}
	else
	{
		const Real supportRadius = particleRadius * static_cast<Real>(4.0);
		float fluidColor[4] = { 0.1f, 0.2f, 0.6f, 1.0f };

		Partio::ParticleAttribute attr;
		if (fluid.partioData->numAttributes() > 0)
			fluid.partioData->attributeInfo(colorField, attr);

		glPointSize(4.0);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < nParticles; i++)
		{
			if (fluid.partioData->numAttributes() > 0)
			{
				float v = 0.0;
				const float* partioVals = NULL;
				if (attr.type == Partio::VECTOR)
				{
					partioVals = fluid.partioData->data<float>(attr, 0);
					v = sqrt(partioVals[3 * i] * partioVals[3 * i] + partioVals[3 * i + 1] * partioVals[3 * i + 1] + partioVals[3 * i + 2] * partioVals[3 * i + 2]);
				}
				else if (attr.type == Partio::FLOAT)
				{
					partioVals = fluid.partioData->data<float>(attr, 0);
					v = partioVals[3 * i];
				}

				v = 0.5f*((v - renderMinValue) / (renderMaxValue - renderMinValue));
				v = min(128.0f*v*v, 0.5f);
				float fluidColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
				MiniGL::hsvToRgb(0.55f, 1.0f, 0.5f + v, fluidColor);
				glColor3fv(fluidColor);
			}
			else
				glColor3fv(fluidColor);
			glVertex3fv(&partioX[3 * i]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}

	float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		pointShaderBegin(&m_shader_vector, particleRadius, &red[0], renderMinValue, renderMaxValue);
		if (fluid.selectedParticles.size() > 0)
		{
			glUniform1f(m_shader_vector.getUniform("radius"), (float)particleRadius*1.05f);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FALSE, GL_FALSE, 0, partioX);

			if (fluid.partioData->numAttributes() > 0)
			{
				Partio::ParticleAttribute attr;
				fluid.partioData->attributeInfo(colorField, attr);
				const float* partioVals = NULL;
				if (attr.type == Partio::VECTOR)
				{
					glEnableVertexAttribArray(1);
					partioVals = fluid.partioData->data<float>(attr, 0);
					glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, &partioVals[0]);
				}
				else if (attr.type == Partio::FLOAT)
				{
					glEnableVertexAttribArray(1);
					partioVals = fluid.partioData->data<float>(attr, 0);
					glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, &partioVals[0]);
				}
			}

			glDrawElements(GL_POINTS, (GLsizei)fluid.selectedParticles.size(), GL_UNSIGNED_INT, fluid.selectedParticles.data());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}
		pointShaderEnd(&m_shader_vector);
	}
	else
	{
		if (fluid.selectedParticles.size() > 0)
		{
			glPointSize(4.0);
			glDisable(GL_LIGHTING);
			glBegin(GL_POINTS);
			for (unsigned int i = 0; i < fluid.selectedParticles.size(); i++)
			{
				glColor3fv(red);
				glVertex3fv(&partioX[3 * fluid.selectedParticles[i]]);
			}
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}
}


void PartioViewer_OpenGL::render(const Boundary &boundary, const bool renderWalls, float *boundaryColor)
{
	if (renderWalls || (!boundary.isWall))
	{
		m_meshShader.begin();
		glUniform1f(m_meshShader.getUniform("shininess"), 5.0f);
		glUniform1f(m_meshShader.getUniform("specular_factor"), 0.2f);

		GLfloat matrix[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
		glUniformMatrix4fv(m_meshShader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
		GLfloat pmatrix[16];
		glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
		glUniformMatrix4fv(m_meshShader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

		glUniform3fv(m_meshShader.getUniform("surface_color"), 1, &boundary.color[0]);

		MiniGL::drawMesh(boundary.mesh, &boundary.color[0]);

		m_meshShader.end();
	}
}

void PartioViewer_OpenGL::hsvToRgb(float h, float s, float v, float *rgb)
{
	int i = (int)floor(h * 6);
	float f = h * 6 - i;
	float p = v * (1 - s);
	float q = v * (1 - f * s);
	float t = v * (1 - (1 - f) * s);

	switch (i % 6)
	{
	case 0: rgb[0] = v, rgb[1] = t, rgb[2] = p; break;
	case 1: rgb[0] = q, rgb[1] = v, rgb[2] = p; break;
	case 2: rgb[0] = p, rgb[1] = v, rgb[2] = t; break;
	case 3: rgb[0] = p, rgb[1] = q, rgb[2] = v; break;
	case 4: rgb[0] = t, rgb[1] = p, rgb[2] = v; break;
	case 5: rgb[0] = v, rgb[1] = p, rgb[2] = q; break;
	}
}

void PartioViewer_OpenGL::updateScalarField()
{
	m_updateScalarField = true;
}

