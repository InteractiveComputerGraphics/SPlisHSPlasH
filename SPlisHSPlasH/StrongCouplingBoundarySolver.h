#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH {
	class StrongCouplingBoundarySolver {
	public:
		~StrongCouplingBoundarySolver(){
			current = nullptr;
		}
		static StrongCouplingBoundarySolver* getCurrent() {
			if (current == nullptr) {
				current = new StrongCouplingBoundarySolver();
			}
			return current;
		}
		void computeDensityAndVolume();
		void computePressureGrad();
		void computeV_s();
		void computeSourceTermRHS();
		void computeSourceTerm();
		void computeDiagonalElement();
	private:
		static StrongCouplingBoundarySolver* current;
	};
}