#ifndef __MiniGL_h__
#define __MiniGL_h__

#include "SPlisHSPlasH/Common.h"
#include <Eigen/Geometry> 
#include "Shader.h"
#include "SPlisHSPlasH/TriangleMesh.h"

#ifdef USE_DOUBLE
#define glNormal3v glNormal3dv
#define glVertex3v glVertex3dv
#define glVertex3 glVertex3d
#define glMultMatrix glMultMatrixd
#define glGetRealv glGetDoublev
#define glLoadMatrix glLoadMatrixd
#define glTranslate glTranslated
#define GL_REAL GL_DOUBLE
#else
#define glNormal3v glNormal3fv
#define glVertex3v glVertex3fv
#define glVertex3 glVertex3f
#define glMultMatrix glMultMatrixf
#define glGetRealv glGetFloatv
#define glLoadMatrix glLoadMatrixf
#define glTranslate glTranslatef
#define GL_REAL GL_FLOAT
#endif


namespace SPH
{
	class MiniGL
	{
	#define IMAGE_ROWS 128
	#define IMAGE_COLS 128
	
	private:
		struct Line
		{
			Vector3r a;
			Vector3r b;
			float color[4];
			float lineWidth;
		};

		struct Point
		{
			Vector3r a;
			float color[4];
			float pointSize;
		};

		struct Triangle
		{
			Vector3r a;
			Vector3r b;
			Vector3r c;
			float color[4];
		};

		struct KeyFunction
		{
			std::function<void()> fct;
			unsigned char key;
		};

		typedef std::function<void()> SceneFct;
		typedef std::function<void()> IdleFct;
		typedef std::function<void(int, int)> ReshapeFct;
		typedef std::function<bool(unsigned char, int, int)> KeyboardFct;
		typedef std::function<bool(int, int, int)> SpecialFct;
		typedef std::function<bool(int, int, int, int)> MousePressFct;
		typedef std::function<bool(int, int)> MouseMoveFct;
		typedef std::function<bool(int, int, int, int)> MouseWheelFct;

		static float fovy;
		static float znear;
		static float zfar;
		static SceneFct scenefunc;
		static void (*exitfunc)(void);
		static IdleFct idlefunc;
		static std::vector<KeyFunction> keyfunc;
		static std::vector<ReshapeFct> m_reshapeFct;
		static std::vector<KeyboardFct> m_keyboardFct;
		static std::vector<SpecialFct> m_specialFct;
		static std::vector<MousePressFct> m_mousePressFct;
		static std::vector<MouseMoveFct> m_mouseMoveFct;
		static std::vector<MouseWheelFct> m_mouseWheelFct;
		static int width;
		static int height;
		static Vector3r m_translation;
		static Quaternionr m_rotation;
		static Real m_zoom;
		static Real movespeed;
		static Real turnspeed;
		static int mouse_button;
		static int modifier_key;
		static int mouse_pos_x_old;
		static int mouse_pos_y_old;
		static int drawMode;
		static unsigned char texData[IMAGE_ROWS][IMAGE_COLS][3];		
		static unsigned int m_texId;
		static void(*selectionfunc) (const Vector2i&, const Vector2i&, void*);
		static void* selectionfuncClientData;
		static void(*mousefunc)(int, int, void*);
		static int mouseFuncButton;		
		static Vector2i m_selectionStart;
		static Real m_quat[4];
		static GLint m_context_major_version;
		static GLint m_context_minor_version;
		static GLint m_context_profile;
		static bool m_breakPointActive;
		static bool m_breakPointLoop;

		static void reshape (int w, int h);
		static void idle ();
		static void keyboard (unsigned char k, int x, int y);
		static void special (int k, int x, int y);
		static void mousePress (int button, int state, int x, int y);
		static void mouseMove (int x, int y);
		static void mouseWheel(int button, int dir, int x, int y);

		static void breakPointMainLoop();
		
	public:
		static void getOpenGLVersion(int &major_version, int &minor_version);
		static void coordinateSystem ();
		static void drawVector(const Vector3r &a, const Vector3r &b, const float w, float *color);
		/** Renders a closed cylinder between two points.
		*/
		static void drawCylinder(const Vector3r &a, const Vector3r &b, const float *color, const float radius = 0.02, const unsigned int subdivisions = 8);
		static void drawSphere(const Vector3r &translation, float radius, float *color, const unsigned int subDivision = 16);
		static void drawTorus(const Vector3r &translation, float innerRadius, float outerRadius, float *color, const unsigned int nsides = 16, const unsigned int rings = 16);
		static void drawQuad (const Vector3r &a, const Vector3r &b, const Vector3r &c, const Vector3r &d, const Vector3r &norm, float *color);
		/** Draw a tetrahedron.
		*/
		static void drawTetrahedron(const Vector3r &a, const Vector3r &b, const Vector3r &c, const Vector3r &d, float *color);
		static void drawTriangle (const Vector3r &a, const Vector3r &b, const Vector3r &c, const Vector3r &norm, float *color);
		static void drawGrid_xz(float *color);
		static void drawGrid_xy(float *color);
		static void drawBitmapText (float x, float y, const char *str, int strLength, float *color);
		static void drawStrokeText(const Real x, const Real y, const Real z, float scale, const char *str, int strLength, float *color);
		static void drawStrokeText (const Vector3r &pos, float scale, const char *str, int strLength, float *color);
		static void drawCube (const Vector3r &translation, const Matrix3r &rotation, float width, float height, float depth, float *color);		
		static void drawPoint (const Vector3r &translation, const float pointSize, const float * const color);
		static void drawMesh(const TriangleMesh &mesh, const float * const color);
		static void drawMesh(const std::vector<Vector3r> &vertices, const std::vector<unsigned int> &faces, const std::vector<Vector3r> &vertexNormals, const float * const color);
		static void setViewport (float pfovy, float pznear, float pzfar, const Vector3r &peyepoint, const Vector3r &plookat);
		static void setViewport (float pfovy, float pznear, float pzfar);
		static void setClientSceneFunc (SceneFct func);
		static void display ();
		static void setClientIdleFunc (IdleFct func);
		static void addKeyFunc(unsigned char k, std::function<void()> func);
		static std::vector<KeyFunction> &getKeyFunc() { return keyfunc; }
		static void init(int argc, char **argv, const int width, const int height, const int posx, const int posy, const char *name);
		static void destroy ();
		static void viewport ();
		static void initLights ();
		static Shader *createShader(const std::string &vertexShader, const std::string &geometryShader, const std::string &fragmentShader);
		static bool checkOpenGLVersion(const int major_version, const int minor_version);
		static void initTexture ();
		static void bindTexture();
		static void unbindTexture();
		static void move (Real x, Real y, Real z);
		static void rotateX (Real x);
		static void rotateY (Real y);
		static void setProjectionMatrix (int width, int height);
		static void setSelectionFunc(void(*func) (const Vector2i&, const Vector2i&, void*), void *clientData);
		static void setMouseMoveFunc(int button, void(*func) (int, int, void*));
		static void unproject(const int x, const int y, Vector3r &pos);
		static float getZNear();
		static float getZFar();
		static void hsvToRgb(float h, float s, float v, float *rgb);
		static void rgbToHsv(float r, float g, float b, float *hsv);
		static int getModifierKey() { return modifier_key; }

		static void addReshapeFunc(ReshapeFct func) { m_reshapeFct.push_back(func); }
		static std::vector<ReshapeFct> &getReshapeFunc() { return m_reshapeFct; }
		static void addKeyboardFunc(KeyboardFct func) { m_keyboardFct.push_back(func); }
		static std::vector<KeyboardFct> &getKeyboardFunc() { return m_keyboardFct; }
		static void addSpecialFunc(SpecialFct func) { m_specialFct.push_back(func); }
		static std::vector<SpecialFct> &getSpecialFunc() { return m_specialFct; }
		static void addMousePressFunc(MousePressFct func) { m_mousePressFct.push_back(func); }
		static std::vector<MousePressFct> &getMousePressFunc() { return m_mousePressFct; }
		static void addMouseMoveFunc(MouseMoveFct func) { m_mouseMoveFct.push_back(func); }
		static std::vector<MouseMoveFct> &getMouseMoveFunc() { return m_mouseMoveFct; }
		static void addMouseWheelFunc(MouseWheelFct func) { m_mouseWheelFct.push_back(func); }
		static std::vector<MouseWheelFct> &getMouseWheelFunc() { return m_mouseWheelFct; }

		static void setBreakPointActive(const bool active);
		static void breakPoint();

		static int getWidth() { return width; }
		static int getHeight() { return height; }
		
		static int getDrawMode() { return drawMode; }
		static void setDrawMode(int val) { drawMode = val; }
		static Quaternionr getRotation() { return m_rotation; }
		static void setRotation(Quaternionr val) { m_rotation = val; }
	};
}

#endif
