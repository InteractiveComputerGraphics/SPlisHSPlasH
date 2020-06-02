#ifndef __PartioViewer_GUI_TweakBar_h__
#define __PartioViewer_GUI_TweakBar_h__

#include "SPlisHSPlasH/Common.h"
#include "../PartioViewer_GUI_Base.h"
#include "extern/AntTweakBar/include/AntTweakBar.h"
#include <vector>

namespace SPH 
{	
	class PartioViewer;

	class PartioViewer_GUI_TweakBar : public PartioViewer_GUI_Base
	{
		public:
			PartioViewer_GUI_TweakBar(PartioViewer *viewer);
			virtual ~PartioViewer_GUI_TweakBar();

		protected:
			PartioViewer *m_viewer;
			unsigned int m_currentFluidModel;
			bool m_renderWalls;
			bool m_showBBox;
			Vector3r m_camPos;
			Vector3r m_camLookat;
			TwBar *m_tweakBar;
			std::map<int, int> m_mapColorField2Attr;

			void initTweakBar();
			void initTweakBarParameters();
			void renderScene();

			static void selection(const Vector2i &start, const Vector2i &end, void *clientData);

			static void TW_CALL setWireframeCB(const void *value, void *clientData);
			static void TW_CALL getWireframeCB(void *value, void *clientData);
			static void TW_CALL setRotationCB(const void *value, void *clientData);
			static void TW_CALL getRotationCB(void *value, void *clientData);

			static void TW_CALL setFrameIndex(const void *value, void *clientData);
			static void TW_CALL getFrameIndex(void *value, void *clientData);
			static void TW_CALL setUsePlane(const void *value, void *clientData);
			static void TW_CALL getUsePlane(void *value, void *clientData);
			static void TW_CALL setPlaneNormal(const void *value, void *clientData);
			static void TW_CALL getPlaneNormal(void *value, void *clientData);
			static void TW_CALL setPlanePoint(const void *value, void *clientData);
			static void TW_CALL getPlanePoint(void *value, void *clientData);
			static void TW_CALL setParticleRadius(const void *value, void *clientData);
			static void TW_CALL getParticleRadius(void *value, void *clientData);
			static void TW_CALL setStartFrame(const void *value, void *clientData);
			static void TW_CALL getStartFrame(void *value, void *clientData);
			static void TW_CALL setEndFrame(const void *value, void *clientData);
			static void TW_CALL getEndFrame(void *value, void *clientData);
			static void TW_CALL setFPS(const void *value, void *clientData);
			static void TW_CALL getFPS(void *value, void *clientData);
			static void TW_CALL setCurrentFluidModel(const void* value, void* clientData);
			static void TW_CALL getCurrentFluidModel(void* value, void* clientData);
					   

		public:
			virtual void init();
			virtual void initParameterGUI();
			virtual void addOptions(cxxopts::Options &options);
			virtual void parseOptions(cxxopts::ParseResult &result);
			virtual void render();
			virtual void update();
			virtual void cleanup();
			virtual void run();
			virtual void stop();
			virtual unsigned int getWidth() const;
			virtual unsigned int getHeight() const;

			TwBar * getTweakBar() const { return m_tweakBar; }
	};
}

#endif