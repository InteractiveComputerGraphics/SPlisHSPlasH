#include "SPlisHSPlasH/Common.h"
#include <Eigen/Dense>
#include <iostream>
#include "GL/glew.h"
#include "Visualization/MiniGL.h"
#include "GL/glut.h"
#include "SPlisHSPlasH/Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/PoissonDiskSampling.h"
#include "Utilities/FileSystem.h"

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

using namespace SPH;
using namespace Eigen;
using namespace std;

void initShader();
void render();
void pointShaderBegin(const float *col);
void pointShaderEnd();
void timeStep() {}


string inputFile = "";
string exePath, dataPath;
Real particleRadius = 0.025;
std::vector<Vector3r> x;
std::vector<Vector3r> v;
Real maxVel = 1.0;
Shader shader;
GLint context_major_version;
GLint context_minor_version;


// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	exePath = FileSystem::getProgramPath();
	dataPath = FileSystem::normalizePath(exePath + "/" + std::string(SPH_DATA_PATH));

	if (argc < 2)
	{
		std::cerr << "Not enough parameters!\n";
		std::cerr << "Usage: PartioViewer.exe [-r radius] particles.bgeo\n";
		return -1;
	}

	bool particleRadiusParam = false;
	for (int i=1; i < argc; i++)
	{
		string argStr = argv[i];
		string type_str = argStr.substr(0, 2);
		if ((type_str == "-r") && (i + 1 < argc))
		{
			particleRadius = stof(argv[++i]);
			particleRadiusParam = true;
		}
		else
			inputFile = argv[i];
	}

	if (!particleRadiusParam)
	{
		particleRadius = 0.025;
		PartioReaderWriter::readParticles(inputFile, Vector3r::Zero(), Matrix3r::Identity(), 1.0, x, v, particleRadius);
	}
	else
		PartioReaderWriter::readParticles(inputFile, Vector3r::Zero(), Matrix3r::Identity(), 1.0, x, v);

	for (unsigned int i = 0; i < v.size(); i++)
		maxVel = std::max(maxVel, v[i].norm());

	// OpenGL
	MiniGL::init(argc, argv, 1024, 768, 0, 0, "Partio Viewer");
	MiniGL::initLights();
	MiniGL::initTweakBarParameters();
	MiniGL::getOpenGLVersion(context_major_version, context_minor_version);
	MiniGL::setViewport(40.0, 0.1f, 500.0, Vector3r(0.0, 3.0, 10.0), Vector3r(0.0, 0.0, 0.0));

	TwAddVarRW(MiniGL::getTweakBar(), "particleRadius", TW_TYPE_REAL, &particleRadius, " label='Particle radius' min=0.001 group=General");

	if (MiniGL::checkOpenGLVersion(3, 3))
		initShader();
	
	MiniGL::setClientSceneFunc(render);
	MiniGL::setClientIdleFunc(1, timeStep);

	glutMainLoop();


	Timing::printAverageTimes();
	Timing::printTimeSums();
	
	return 0;
}


void initShader()
{
	string vertFile = dataPath + "/shaders/vs_points.glsl";
	string fragFile = dataPath + "/shaders/fs_points.glsl";
	shader.compileShaderFile(GL_VERTEX_SHADER, vertFile);
	shader.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
	shader.createAndLinkProgram();
	shader.begin();
	shader.addUniform("modelview_matrix");
	shader.addUniform("projection_matrix");
	shader.addUniform("radius");
	shader.addUniform("viewport_width");
	shader.addUniform("color");
	shader.addUniform("projection_radius");
	shader.addUniform("max_velocity");
	shader.end();
}

void render()
{
	MiniGL::coordinateSystem();

	// Draw simulation model
	const unsigned int nParticles = (unsigned int) x.size();

	if (MiniGL::checkOpenGLVersion(3, 3))
	{		
		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		pointShaderBegin(&fluidColor[0]);

		if (nParticles > 0)
		{
			glUniform1f(shader.getUniform("max_velocity"), (GLfloat) maxVel);

			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, x.data());
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, v.data());
			glDrawArrays(GL_POINTS, 0, nParticles);
			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}

		pointShaderEnd();
	}
	else
	{
		const Real supportRadius = particleRadius*4.0;
		float fluidColor[4] = { 0.1f, 0.2f, 0.6f, 1.0f };

		glPointSize(4.0);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < nParticles; i++)
		{
			// modify color according to the velocity
			Eigen::Vector3f hsv;
			MiniGL::hsvToRgb(fluidColor[0], fluidColor[1], fluidColor[2], &hsv[0]);
			Real vl = v[i].norm();
			vl = std::min((1.0 / maxVel)*vl, 1.0);
			float finalColor[3];
			MiniGL::hsvToRgb(hsv[0], std::max(1.0f - (float) vl, 0.0f), 1.0f, finalColor);
			glColor3fv(finalColor);
			glVertex3v(&x[i][0]);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}
}

void pointShaderBegin(const float *col)
{
	shader.begin();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glUniform1f(shader.getUniform("viewport_width"), (float)viewport[2]);
	glUniform1f(shader.getUniform("radius"), (float)particleRadius);
	glUniform3fv(shader.getUniform("color"), 1, col);

	GLfloat matrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
	glUniformMatrix4fv(shader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
	GLfloat pmatrix[16];
	glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
	glUniformMatrix4fv(shader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

	glEnable(GL_DEPTH_TEST);
	// Point sprites do not have to be explicitly enabled since OpenGL 3.2 where
	// they are enabled by default. Moreover GL_POINT_SPRITE is deprecate and only
	// supported before OpenGL 3.2 or with compatibility profile enabled.
	glEnable(GL_POINT_SPRITE);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glPointParameterf(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);
}

void pointShaderEnd()
{
	shader.end();
}
