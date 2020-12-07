#include "Simulator_OpenGL.h"

#include "GUI/OpenGL/MiniGL.h"
#include "GUI/OpenGL/colormaps/colormap_jet.h"
#include "GUI/OpenGL/colormaps/colormap_plasma.h"
#include "GUI/OpenGL/colormaps/colormap_bwr.h"
#include "GUI/OpenGL/colormaps/colormap_coolwarm.h"
#include "GUI/OpenGL/colormaps/colormap_seismic.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include <Utilities/Timing.h>

using namespace SPH;
using namespace std;

Shader Simulator_OpenGL::m_shader_scalar;
Shader Simulator_OpenGL::m_shader_scalar_map;
Shader Simulator_OpenGL::m_shader_vector;
Shader Simulator_OpenGL::m_meshShader;
GLuint Simulator_OpenGL::m_textureMap = 0;


Simulator_OpenGL::Simulator_OpenGL()
{	
}

Simulator_OpenGL::~Simulator_OpenGL(void)
{	
}


void Simulator_OpenGL::initShaders(const std::string &shaderPath)
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


void Simulator_OpenGL::meshShaderBegin(const float *col)
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
	glUniform3fv(m_meshShader.getUniform("surface_color"), 1, col);
}

void Simulator_OpenGL::meshShaderEnd()
{
	m_meshShader.end();
}

void Simulator_OpenGL::pointShaderBegin(Shader *shader, const Real particleRadius, const float *col, const Real minVal, const Real maxVal, const bool useTexture, float const* color_map)
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


void Simulator_OpenGL::pointShaderEnd(Shader *shader, const bool useTexture)
{
	glBindTexture(GL_TEXTURE_1D, 0);
	shader->end();
}

void Simulator_OpenGL::renderFluid(FluidModel *model, float *fluidColor,
	const unsigned int colorMapType, const bool useScalarField, const std::vector<float> &scalarField,
	const Real renderMinValue, const Real renderMaxValue)
{
	// Draw simulation model
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nParticles = model->numActiveParticles();
	if (nParticles == 0)
		return;

	float surfaceColor[4] = { 0.2f, 0.6f, 0.8f, 1 };
	float speccolor[4] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
	glColor3fv(surfaceColor);

	const Real particleRadius = sim->getParticleRadius();
	const Real supportRadius = sim->getSupportRadius();
	Real vmax = static_cast<Real>(0.4*2.0)*supportRadius / TimeManager::getCurrent()->getTimeStepSize();
	Real vmin = 0.0;

	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		Shader *shader_scalar = &m_shader_scalar_map;
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
			shader_scalar = &m_shader_scalar;

		if (!useScalarField)
			pointShaderBegin(shader_scalar, particleRadius, &fluidColor[0], renderMinValue, renderMaxValue, false);
		else 
			pointShaderBegin(shader_scalar, particleRadius, &fluidColor[0], renderMinValue, renderMaxValue, true, color_map);

		if (model->numActiveParticles() > 0)
		{
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &model->getPosition(0));

			if (useScalarField)
			{
				glEnableVertexAttribArray(1);
				glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, &scalarField[0]);
			}

			glDrawArrays(GL_POINTS, 0, model->numActiveParticles());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}

		if (!useScalarField)
			pointShaderEnd(shader_scalar, false);
		else 
			pointShaderEnd(shader_scalar, true);
	}
	else
	{
		glPointSize(4.0);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < nParticles; i++)
		{
			Real v = model->getVelocity(i).norm();
			v = static_cast<Real>(0.5)*((v - vmin) / (vmax - vmin));
			v = min(static_cast<Real>(128.0)*v*v, static_cast<Real>(0.5));
			float fluidColor[4] = { 0.2f, 0.2f, 0.2f, 1.0 };
			MiniGL::hsvToRgb(0.55f, 1.0f, 0.5f + (float)v, fluidColor);

			glColor3fv(fluidColor);
			glVertex3v(&model->getPosition(i)[0]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}
}

void Simulator_OpenGL::renderSelectedParticles(FluidModel *model, const std::vector<std::vector<unsigned int>>& selectedParticles,
		const unsigned int colorMapType, 
		const Real renderMinValue, const Real renderMaxValue)
{
	float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
	const unsigned int fluidIndex = model->getPointSetIndex();
	Simulation *sim = Simulation::getCurrent();
	const Real particleRadius = sim->getParticleRadius();
	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		if ((selectedParticles.size() > 0) && ((selectedParticles[fluidIndex].size() > 0)))
		{
			pointShaderBegin(&m_shader_vector, particleRadius, &red[0], renderMinValue, renderMaxValue);
			const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
			glUniform1f(m_shader_vector.getUniform("radius"), (float)radius*1.05f);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &model->getPosition(0));
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_REAL, GL_FALSE, 0, &model->getVelocity(0));
			glDrawElements(GL_POINTS, (GLsizei)selectedParticles[fluidIndex].size(), GL_UNSIGNED_INT, selectedParticles[fluidIndex].data());
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
			pointShaderEnd(&m_shader_vector);
		}
	}
	else
	{
		if (selectedParticles.size() > 0)
		{
			glPointSize(4.0);
			glDisable(GL_LIGHTING);
			glBegin(GL_POINTS);
			for (unsigned int i = 0; i < selectedParticles[fluidIndex].size(); i++)
			{
				glColor3fv(red);
				glVertex3v(&model->getPosition(selectedParticles[fluidIndex][i])[0]);
			}
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}

}

void Simulator_OpenGL::renderBoundaryParticles(BoundaryModel_Akinci2012 *model, const float *col)
{
	Simulation *sim = Simulation::getCurrent();
	const Real particleRadius = sim->getParticleRadius();

	if (MiniGL::checkOpenGLVersion(3, 3))
	{
		Simulator_OpenGL::pointShaderBegin(&m_shader_scalar, particleRadius, col, 0.0, 100000.0);
		glEnableVertexAttribArray(0);

		glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &model->getPosition(0)[0]);
		glDrawArrays(GL_POINTS, 0, model->numberOfParticles());
		glDisableVertexAttribArray(0);

		Simulator_OpenGL::pointShaderEnd(&m_shader_scalar);
	}
	else
	{
		glDisable(GL_LIGHTING);
		glPointSize(4.0f);

		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < model->numberOfParticles(); i++)
		{
			glColor3fv(col);
			glVertex3v(&model->getPosition(i)[0]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}
}

void Simulator_OpenGL::renderBoundary(BoundaryModel *model, const float *col)
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

	glUniform3fv(m_meshShader.getUniform("surface_color"), 1, col);

	const std::vector<Vector3r> &vertices = model->getRigidBodyObject()->getVertices();
	const std::vector<Vector3r> &vNormals = model->getRigidBodyObject()->getVertexNormals();
	const std::vector<unsigned int> &faces = model->getRigidBodyObject()->getFaces();

	MiniGL::drawMesh(vertices, faces, vNormals, col);

	m_meshShader.end();
}