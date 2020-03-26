#ifndef __Simulator_GUI_TweakBar_h__
#define __Simulator_GUI_TweakBar_h__

#include "SPlisHSPlasH/Common.h"
#include "../Simulator_GUI_Base.h"
#include "extern/AntTweakBar/include/AntTweakBar.h"
#include <vector>

namespace SPH 
{	
	class Simulator_GUI_TweakBar : public Simulator_GUI_Base
	{
		public:
			Simulator_GUI_TweakBar(SimulatorBase *base);
			virtual ~Simulator_GUI_TweakBar();

		protected:
			unsigned int m_currentFluidModel;
			TwBar *m_tweakBar;
			Vector3r m_oldMousePos;
			std::vector<std::string> m_colorFieldNames;
			std::vector<std::vector<unsigned int>> m_selectedParticles;			

			std::vector<std::vector<unsigned int>>& getSelectedParticles() { return m_selectedParticles; }
			void initTweakBar();
			void initTweakBarParameters();
			void renderBoundary();
			void particleInfo();

			static void selection(const Vector2i &start, const Vector2i &end, void *clientData);
			static void mouseMove(int x, int y, void *clientData);

			static void TW_CALL setTimeCB(const void *value, void *clientData);
			static void TW_CALL getTimeCB(void *value, void *clientData);
			static void TW_CALL setWireframeCB(const void *value, void *clientData);
			static void TW_CALL getWireframeCB(void *value, void *clientData);
			static void TW_CALL setRotationCB(const void *value, void *clientData);
			static void TW_CALL getRotationCB(void *value, void *clientData);
			static void TW_CALL setCurrentFluidModel(const void *value, void *clientData);
			static void TW_CALL getCurrentFluidModel(void *value, void *clientData);
			static void TW_CALL setColorField(const void *value, void *clientData);
			static void TW_CALL getColorField(void *value, void *clientData);
			static void TW_CALL setRenderMaxValue(const void *value, void *clientData);
			static void TW_CALL getRenderMaxValue(void *value, void *clientData);
			static void TW_CALL setRenderMinValue(const void *value, void *clientData);
			static void TW_CALL getRenderMinValue(void *value, void *clientData);
			static void TW_CALL setColorMapType(const void *value, void *clientData);
			static void TW_CALL getColorMapType(void *value, void *clientData);

		public:
			virtual void init(int argc, char **argv, const char *name);
			virtual void initParameterGUI();
			virtual void initSimulationParameterGUI();
			virtual void render();
			virtual void reset();
			virtual void update();
			virtual void cleanup();
			virtual void run();
			virtual void stop();
			virtual void addKeyFunc(char k, std::function<void()> const& func);

			TwBar * getTweakBar() const { return m_tweakBar; }
	};
}

#endif