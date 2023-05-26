#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH {
	class StrongCouplingBoundarySolver {
	private:
		// values required for Gissler 2019 strong coupling based on Akinci 2012
		std::vector<std::vector<double>> m_density;
		std::vector<std::vector<double>> m_pressure;
		std::vector<std::vector<double>> m_lastPressure;
		std::vector<std::vector<double>> m_artificialVolume;
		std::vector<std::vector<Vector3r>> m_v_s;
		std::vector<std::vector<double>> m_s; // source term
		std::vector<std::vector<Vector3r>> m_pressureGrad;
		std::vector<std::vector<Vector3r>> m_v_rr;
		std::vector<std::vector<Vector3r>> m_v_b;
		std::vector<std::vector<Vector3r>> m_grad_p_b;
		std::vector<std::vector<Vector3r>> m_predictVelocity;
		std::vector<std::vector<Vector3r>> m_predictPosition;
		std::vector<std::vector<double>> m_minus_rho_div_v_rr; // RHS to the source term
		std::vector<std::vector<double>> m_diagonalElement; // diagonal element for jacobi iteration
		std::vector<double> m_restDensity;
		std::vector<Vector3r> m_v_rr_body;
		std::vector<Vector3r> m_omega_rr_body;
		std::vector<Vector3r> m_v_b_body;
		std::vector<Vector3r> m_omega_b_body;

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
		void reset();
		void computeDensityAndVolume();
		void computePressureGrad();
		void computeV_s();
		void computeSourceTermRHS();
		void computeSourceTerm();
		void computeDiagonalElement();
		void computeFriction();
		void performNeighborhoodSearchSort();

		FORCE_INLINE const double& getDensity(const unsigned int& rigidBodyIndex, const unsigned int& index) const {
			return m_density[rigidBodyIndex][index];
		}

		FORCE_INLINE double& getDensity(const unsigned int& rigidBodyIndex, const unsigned int& index) {
			return m_density[rigidBodyIndex][index];
		}

		FORCE_INLINE void setDensity(const unsigned int& rigidBodyIndex, const unsigned int& index, const double& value) {
			m_density[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const double& getPressure(const unsigned int& rigidBodyIndex, const unsigned int& index) const {
			return m_pressure[rigidBodyIndex][index];
		}

		FORCE_INLINE double& getPressure(const unsigned int& rigidBodyIndex, const unsigned int& index) {
			return m_pressure[rigidBodyIndex][index];
		}

		FORCE_INLINE void setPressure(const unsigned int& rigidBodyIndex, const unsigned int& index, const double& value) {
			m_pressure[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const double& getLastPressure(const unsigned int& rigidBodyIndex, const unsigned int& index) const {
			return m_lastPressure[rigidBodyIndex][index];
		}

		FORCE_INLINE double& getLastPressure(const unsigned int& rigidBodyIndex, const unsigned int& index) {
			return m_lastPressure[rigidBodyIndex][index];
		}

		FORCE_INLINE void setLastPressure(const unsigned int& rigidBodyIndex, const unsigned int& index, const double& value) {
			m_lastPressure[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const double& getRestDensity(const unsigned int& rigidBodyIndex) const {
			return m_restDensity[rigidBodyIndex];
		}

		FORCE_INLINE double& getRestDensity(const unsigned int& rigidBodyIndex) {
			return m_restDensity[rigidBodyIndex];
		}

		FORCE_INLINE void setRestDensity(const unsigned int& rigidBodyIndex, const double& value) {
			m_restDensity[rigidBodyIndex] = value;
		}

		FORCE_INLINE const Vector3r& getV_rr_body(const unsigned int& rigidBodyIndex) const {
			return m_v_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE Vector3r& getV_rr_body(const unsigned int& rigidBodyIndex) {
			return m_v_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE void setV_rr_body(const unsigned int& rigidBodyIndex, const Vector3r& value) {
			m_v_rr_body[rigidBodyIndex] = value;
		}

		FORCE_INLINE const Vector3r& getV_b_body(const unsigned int& rigidBodyIndex) const {
			return m_v_b_body[rigidBodyIndex];
		}

		FORCE_INLINE Vector3r& getV_b_body(const unsigned int& rigidBodyIndex) {
			return m_v_b_body[rigidBodyIndex];
		}

		FORCE_INLINE void setV_b_body(const unsigned int& rigidBodyIndex, const Vector3r& value) {
			m_v_b_body[rigidBodyIndex] = value;
		}

		FORCE_INLINE const Vector3r& getOmega_rr_body(const unsigned int& rigidBodyIndex) const {
			return m_omega_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE Vector3r& getOmega_rr_body(const unsigned int& rigidBodyIndex) {
			return m_omega_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE void setOmega_rr_body(const unsigned int& rigidBodyIndex, const Vector3r& value) {
			m_omega_rr_body[rigidBodyIndex] = value;
		}

		FORCE_INLINE const Vector3r& getOmega_b_body(const unsigned int& rigidBodyIndex) const {
			return m_omega_b_body[rigidBodyIndex];
		}

		FORCE_INLINE Vector3r& getOmega_b_body(const unsigned int& rigidBodyIndex) {
			return m_omega_b_body[rigidBodyIndex];
		}

		FORCE_INLINE void setOmega_b_body(const unsigned int& rigidBodyIndex, const Vector3r& value) {
			m_omega_b_body[rigidBodyIndex] = value;
		}
		
		FORCE_INLINE Vector3r& getV_s(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_v_s[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getV_s(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_v_s[rigidBodyIndex][i];
		}

		FORCE_INLINE void setV_s(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_v_s[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getV_rr(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getV_rr(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE void setV_rr(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_v_rr[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getV_b(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_v_b[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getV_b(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_v_b[rigidBodyIndex][i];
		}

		FORCE_INLINE void setV_b(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_v_b[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getGrad_p_b(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_grad_p_b[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getGrad_p_b(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_grad_p_b[rigidBodyIndex][i];
		}

		FORCE_INLINE void setGrad_p_b(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_grad_p_b[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getPredictVelocity(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_predictVelocity[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getPredictVelocity(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_predictVelocity[rigidBodyIndex][i];
		}

		FORCE_INLINE void setPredictVelocity(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_predictVelocity[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getPredictPosition(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_predictPosition[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getPredictPosition(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_predictPosition[rigidBodyIndex][i];
		}

		FORCE_INLINE void setPredictPosition(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_predictPosition[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getPressureGrad(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_pressureGrad[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getPressureGrad(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_pressureGrad[rigidBodyIndex][i];
		}

		FORCE_INLINE void setPressureGrad(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_pressureGrad[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE double& getSourceTermRHS(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_minus_rho_div_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE const double& getSourceTermRHS(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_minus_rho_div_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE void setSourceTermRHS(const unsigned int& rigidBodyIndex, const unsigned int i, const double& value) {
			m_minus_rho_div_v_rr[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE double& getDiagonalElement(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_diagonalElement[rigidBodyIndex][i];
		}

		FORCE_INLINE const double& getDiagonalElement(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_diagonalElement[rigidBodyIndex][i];
		}

		FORCE_INLINE void setDiagonalElement(const unsigned int& rigidBodyIndex, const unsigned int i, const double& value) {
			m_diagonalElement[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE double& getSourceTerm(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_s[rigidBodyIndex][i];
		}

		FORCE_INLINE const double& getSourceTerm(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_s[rigidBodyIndex][i];
		}

		FORCE_INLINE void setSourceTerm(const unsigned int& rigidBodyIndex, const unsigned int i, const double& value) {
			m_s[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE const double& getArtificialVolume(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_artificialVolume[rigidBodyIndex][i];
		}

		FORCE_INLINE double& getArtificialVolume(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_artificialVolume[rigidBodyIndex][i];
		}

		FORCE_INLINE void setArtificialVolume(const unsigned int& rigidBodyIndex, const unsigned int i, const double& val) {
			m_artificialVolume[rigidBodyIndex][i] = val;
		}
	};
}