#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH {
	class StrongCouplingBoundarySolver {
	private:
		// values required for Gissler 2019 strong coupling based on Akinci 2012
		std::vector<std::vector<Real>> m_density;
		std::vector<std::vector<Real>> m_pressure;
		std::vector<std::vector<Real>> m_lastPressure;
		std::vector<std::vector<Real>> m_artificialVolume;
		std::vector<std::vector<Vector3r>> m_v_s;
		std::vector<std::vector<Real>> m_s; // source term
		std::vector<std::vector<Vector3r>> m_pressureGrad;
		std::vector<std::vector<Vector3r>> m_v_rr;
		std::vector<std::vector<Vector3r>> m_predictVelocity;
		std::vector<std::vector<Real>> m_minus_rho_div_v_rr; // RHS to the source term
		std::vector<std::vector<Real>> m_diagonalElement; // diagonal element for jacobi iteration
		std::vector<Real> m_restDensity;
		std::vector<Vector3r> m_v_rr_body;
		std::vector<Vector3r> m_omega_rr_body;

		static StrongCouplingBoundarySolver* current;

	public:
		StrongCouplingBoundarySolver();
		~StrongCouplingBoundarySolver();
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
		void performNeighborhoodSearchSort();

		FORCE_INLINE const Real& getDensity(const Real& rigidBodyIndex, const Real& index) const {
			return m_density[rigidBodyIndex][index];
		}

		FORCE_INLINE Real& getDensity(const Real& rigidBodyIndex, const Real& index) {
			return m_density[rigidBodyIndex][index];
		}

		FORCE_INLINE void setDensity(const Real& rigidBodyIndex, const Real& index, const Real& value) {
			m_density[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const Real& getPressure(const Real& rigidBodyIndex, const Real& index) const {
			return m_pressure[rigidBodyIndex][index];
		}

		FORCE_INLINE Real& getPressure(const Real& rigidBodyIndex, const Real& index) {
			return m_pressure[rigidBodyIndex][index];
		}

		FORCE_INLINE void setPressure(const Real& rigidBodyIndex, const Real& index, const Real& value) {
			m_pressure[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const Real& getLastPressure(const Real& rigidBodyIndex, const Real& index) const {
			return m_lastPressure[rigidBodyIndex][index];
		}

		FORCE_INLINE Real& getLastPressure(const Real& rigidBodyIndex, const Real& index) {
			return m_lastPressure[rigidBodyIndex][index];
		}

		FORCE_INLINE void setLastPressure(const Real& rigidBodyIndex, const Real& index, const Real& value) {
			m_lastPressure[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const Real& getRestDensity(const Real& rigidBodyIndex) const {
			return m_restDensity[rigidBodyIndex];
		}

		FORCE_INLINE Real& getRestDensity(const Real& rigidBodyIndex) {
			return m_restDensity[rigidBodyIndex];
		}

		FORCE_INLINE void setRestDensity(const Real& rigidBodyIndex, const Real& value) {
			m_restDensity[rigidBodyIndex] = value;
		}

		FORCE_INLINE const Vector3r& getV_rr_body(const Real& rigidBodyIndex) const {
			return m_v_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE Vector3r& getV_rr_body(const Real& rigidBodyIndex) {
			return m_v_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE void setV_rr_body(const Real& rigidBodyIndex, const Vector3r& value) {
			m_v_rr_body[rigidBodyIndex] = value;
		}

		FORCE_INLINE const Vector3r& getOmega_rr_body(const Real& rigidBodyIndex) const {
			return m_omega_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE Vector3r& getOmega_rr_body(const Real& rigidBodyIndex) {
			return m_omega_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE void setOmega_rr_body(const Real& rigidBodyIndex, const Vector3r& value) {
			m_omega_rr_body[rigidBodyIndex] = value;
		}
		
		FORCE_INLINE Vector3r& getV_s(const Real& rigidBodyIndex, const unsigned int i) {
			return m_v_s[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getV_s(const Real& rigidBodyIndex, const unsigned int i) const {
			return m_v_s[rigidBodyIndex][i];
		}

		FORCE_INLINE void setV_s(const Real& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_v_s[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getV_rr(const Real& rigidBodyIndex, const unsigned int i) {
			return m_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getV_rr(const Real& rigidBodyIndex, const unsigned int i) const {
			return m_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE void setV_rr(const Real& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_v_rr[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getPredictVelocity(const Real& rigidBodyIndex, const unsigned int i) {
			return m_predictVelocity[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getPredictVelocity(const Real& rigidBodyIndex, const unsigned int i) const {
			return m_predictVelocity[rigidBodyIndex][i];
		}

		FORCE_INLINE void setPredictVelocity(const Real& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_predictVelocity[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getPressureGrad(const Real& rigidBodyIndex, const unsigned int i) {
			return m_pressureGrad[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getPressureGrad(const Real& rigidBodyIndex, const unsigned int i) const {
			return m_pressureGrad[rigidBodyIndex][i];
		}

		FORCE_INLINE void setPressureGrad(const Real& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_pressureGrad[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Real& getSourceTermRHS(const Real& rigidBodyIndex, const unsigned int i) {
			return m_minus_rho_div_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE const Real& getSourceTermRHS(const Real& rigidBodyIndex, const unsigned int i) const {
			return m_minus_rho_div_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE void setSourceTermRHS(const Real& rigidBodyIndex, const unsigned int i, const Real& value) {
			m_minus_rho_div_v_rr[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Real& getDiagonalElement(const Real& rigidBodyIndex, const unsigned int i) {
			return m_diagonalElement[rigidBodyIndex][i];
		}

		FORCE_INLINE const Real& getDiagonalElement(const Real& rigidBodyIndex, const unsigned int i) const {
			return m_diagonalElement[rigidBodyIndex][i];
		}

		FORCE_INLINE void setDiagonalElement(const Real& rigidBodyIndex, const unsigned int i, const Real& value) {
			m_diagonalElement[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Real& getSourceTerm(const Real& rigidBodyIndex, const unsigned int i) {
			return m_s[rigidBodyIndex][i];
		}

		FORCE_INLINE const Real& getSourceTerm(const Real& rigidBodyIndex, const unsigned int i) const {
			return m_s[rigidBodyIndex][i];
		}

		FORCE_INLINE void setSourceTerm(const Real& rigidBodyIndex, const unsigned int i, const Real& value) {
			m_s[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE const Real& getArtificialVolume(const Real& rigidBodyIndex, const unsigned int i) const {
			return m_artificialVolume[rigidBodyIndex][i];
		}

		FORCE_INLINE Real& getArtificialVolume(const Real& rigidBodyIndex, const unsigned int i) {
			return m_artificialVolume[rigidBodyIndex][i];
		}

		FORCE_INLINE void setArtificialVolume(const Real& rigidBodyIndex, const unsigned int i, const Real& val) {
			m_artificialVolume[rigidBodyIndex][i] = val;
		}
	};
}