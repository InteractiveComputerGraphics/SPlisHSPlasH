#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "ParameterObject.h"
#include <algorithm>

namespace SPH {
	/** \brief This class stores the information of a dynamic rigid body which
	* is used for the strong coupling method introduced in
	* Interlinked SPH Pressure Solvers for Strong Fluid-Rigid Coupling. Gissler et al. https://doi.org/10.1145/3284980
	*/
	class StrongCouplingBoundarySolver : public GenParam::ParameterObject {
	private:
		// values required for Gissler 2019 strong coupling based on Akinci 2012
		Real m_restDensity;
		std::vector<std::vector<Real>> m_density;
		std::vector<std::vector<Real>> m_pressure;
		std::vector<std::vector<Real>> m_lastPressure;
		std::vector<std::vector<Real>> m_artificialVolume;
		std::vector<std::vector<Vector3r>> m_v_s;
		std::vector<std::vector<Real>> m_s; // source term
		std::vector<std::vector<Vector3r>> m_pressureGrad;
		std::vector<std::vector<Vector3r>> m_v_rr;
		std::vector<std::vector<Vector3r>> m_predictedVelocity;
		std::vector<std::vector<Vector3r>> m_predictedPosition;
		std::vector<std::vector<Real>> m_minus_rho_div_v_s;
		std::vector<std::vector<Real>> m_minus_rho_div_v_rr; // RHS to the source term
		std::vector<std::vector<Real>> m_diagonalElement; // diagonal element for jacobi iteration
		std::vector<Vector3r> m_v_rr_body;
		std::vector<Vector3r> m_omega_rr_body;
		std::vector<std::vector<unsigned int>> m_contacts;
		std::vector<unsigned int> m_bodyContacts;
		unsigned int m_contactsAllBodies;

		static StrongCouplingBoundarySolver* current;

		// only used in WCSPH
		unsigned int m_maxIterations;
		unsigned int m_minIterations;

		Real m_maxDensityDeviation;

		int m_kernelMethod;
		int m_gradKernelMethod;
		Real(*m_kernelFct)(const Vector3r&);
		Vector3r(*m_gradKernelFct)(const Vector3r& r);
		Real m_W_zero;

	public:
		static int KERNEL_METHOD;
		static int GRAD_KERNEL_METHOD;
		static int ENUM_KERNEL_CUBIC;
		static int ENUM_KERNEL_WENDLANDQUINTICC2;
		static int ENUM_KERNEL_POLY6;
		static int ENUM_KERNEL_SPIKY;
		static int ENUM_KERNEL_CUBIC_2D;
		static int ENUM_KERNEL_WENDLANDQUINTICC2_2D;
		static int ENUM_GRADKERNEL_CUBIC;
		static int ENUM_GRADKERNEL_WENDLANDQUINTICC2;
		static int ENUM_GRADKERNEL_POLY6;
		static int ENUM_GRADKERNEL_SPIKY;
		static int ENUM_GRADKERNEL_CUBIC_2D;
		static int ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D;

		StrongCouplingBoundarySolver();
		~StrongCouplingBoundarySolver();
		static StrongCouplingBoundarySolver* getCurrent();
		virtual void initParameters();
		void reset();
		void resize(unsigned int size);
		void computeContacts();
		void computeDensityAndVolume();
		void computeV_s();
		void computeSourceTermRHS();
		void computeSourceTermRHSForBody(const unsigned int& boundaryPointSetIndex);
		void computeSourceTerm();
		void computeDiagonalElement();
		void pressureSolveIteration(Real& avgDensityDeviation);
		Vector3r copmuteParticleFriction(const unsigned int& boundaryPointSetIndex, const unsigned int& index, const Vector3r& force);
		void applyForce();
		void performNeighborhoodSearchSort();
		Real W(const Vector3r& r);
		Vector3r gradW(const Vector3r& r);
		Real W_zero();

		

		FORCE_INLINE const Real& getDensity(const unsigned int& rigidBodyIndex, const unsigned int& index) const {
			return m_density[rigidBodyIndex][index];
		}

		FORCE_INLINE Real& getDensity(const unsigned int& rigidBodyIndex, const unsigned int& index) {
			return m_density[rigidBodyIndex][index];
		}

		FORCE_INLINE void setDensity(const unsigned int& rigidBodyIndex, const unsigned int& index, const Real& value) {
			m_density[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const Real& getPressure(const unsigned int& rigidBodyIndex, const unsigned int& index) const {
			return m_pressure[rigidBodyIndex][index];
		}

		FORCE_INLINE Real& getPressure(const unsigned int& rigidBodyIndex, const unsigned int& index) {
			return m_pressure[rigidBodyIndex][index];
		}

		FORCE_INLINE void setPressure(const unsigned int& rigidBodyIndex, const unsigned int& index, const Real& value) {
			m_pressure[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const Real& getLastPressure(const unsigned int& rigidBodyIndex, const unsigned int& index) const {
			return m_lastPressure[rigidBodyIndex][index];
		}

		FORCE_INLINE Real& getLastPressure(const unsigned int& rigidBodyIndex, const unsigned int& index) {
			return m_lastPressure[rigidBodyIndex][index];
		}

		FORCE_INLINE void setLastPressure(const unsigned int& rigidBodyIndex, const unsigned int& index, const Real& value) {
			m_lastPressure[rigidBodyIndex][index] = value;
		}

		FORCE_INLINE const Real& getRestDensity() const {
			return m_restDensity;
		}

		FORCE_INLINE Real& getRestDensity() {
			return m_restDensity;
		}

		FORCE_INLINE void setRestDensity(const Real& value) {
			m_restDensity = value;
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

		FORCE_INLINE const Vector3r& getOmega_rr_body(const unsigned int& rigidBodyIndex) const {
			return m_omega_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE Vector3r& getOmega_rr_body(const unsigned int& rigidBodyIndex) {
			return m_omega_rr_body[rigidBodyIndex];
		}

		FORCE_INLINE void setOmega_rr_body(const unsigned int& rigidBodyIndex, const Vector3r& value) {
			m_omega_rr_body[rigidBodyIndex] = value;
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

		FORCE_INLINE Vector3r& getPredictedVelocity(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_predictedVelocity[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getPredictedVelocity(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_predictedVelocity[rigidBodyIndex][i];
		}

		FORCE_INLINE void setPredictedVelocity(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_predictedVelocity[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Vector3r& getPredictedPosition(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_predictedPosition[rigidBodyIndex][i];
		}

		FORCE_INLINE const Vector3r& getPredictedPosition(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_predictedPosition[rigidBodyIndex][i];
		}

		FORCE_INLINE void setPredictedPosition(const unsigned int& rigidBodyIndex, const unsigned int i, const Vector3r& value) {
			m_predictedPosition[rigidBodyIndex][i] = value;
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

		FORCE_INLINE Real& getMinus_rho_div_v_s(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_minus_rho_div_v_s[rigidBodyIndex][i];
		}

		FORCE_INLINE const Real& getMinus_rho_div_v_s(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_minus_rho_div_v_s[rigidBodyIndex][i];
		}

		FORCE_INLINE void setMinus_rho_div_v_s(const unsigned int& rigidBodyIndex, const unsigned int i, const Real& value) {
			m_minus_rho_div_v_s[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Real& getSourceTermRHS(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_minus_rho_div_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE const Real& getSourceTermRHS(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_minus_rho_div_v_rr[rigidBodyIndex][i];
		}

		FORCE_INLINE void setSourceTermRHS(const unsigned int& rigidBodyIndex, const unsigned int i, const Real& value) {
			m_minus_rho_div_v_rr[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Real& getDiagonalElement(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_diagonalElement[rigidBodyIndex][i];
		}

		FORCE_INLINE const Real& getDiagonalElement(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_diagonalElement[rigidBodyIndex][i];
		}

		FORCE_INLINE void setDiagonalElement(const unsigned int& rigidBodyIndex, const unsigned int i, const Real& value) {
			m_diagonalElement[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE Real& getSourceTerm(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_s[rigidBodyIndex][i];
		}

		FORCE_INLINE const Real& getSourceTerm(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_s[rigidBodyIndex][i];
		}

		FORCE_INLINE void setSourceTerm(const unsigned int& rigidBodyIndex, const unsigned int i, const Real& value) {
			m_s[rigidBodyIndex][i] = value;
		}

		FORCE_INLINE const Real& getArtificialVolume(const unsigned int& rigidBodyIndex, const unsigned int i) const {
			return m_artificialVolume[rigidBodyIndex][i];
		}

		FORCE_INLINE Real& getArtificialVolume(const unsigned int& rigidBodyIndex, const unsigned int i) {
			return m_artificialVolume[rigidBodyIndex][i];
		}

		FORCE_INLINE void setArtificialVolume(const unsigned int& rigidBodyIndex, const unsigned int i, const Real& val) {
			m_artificialVolume[rigidBodyIndex][i] = val;
		}

		FORCE_INLINE const unsigned int& getBodyContacts(const unsigned int& rigidBodyIndex) const {
			return m_bodyContacts[rigidBodyIndex];
		}

		FORCE_INLINE unsigned int& getBodyContacts(const unsigned int& rigidBodyIndex) {
			return m_bodyContacts[rigidBodyIndex];
		}

		FORCE_INLINE void setBodyContacts(const unsigned int& rigidBodyIndex, const unsigned int& value) {
			m_bodyContacts[rigidBodyIndex] = value;
		}

		FORCE_INLINE unsigned int getAllContacts() {
			return m_contactsAllBodies;
		}

		FORCE_INLINE const unsigned int& getMaxIterations() const {
			return m_maxIterations;
		}
		FORCE_INLINE unsigned int& getMaxIterations() {
			return m_maxIterations;
		}
		FORCE_INLINE void setMaxIterations(const unsigned int& value) {
			m_maxIterations = value;
		}

		FORCE_INLINE const unsigned int& getMinIterations() const {
			return m_minIterations;
		}
		FORCE_INLINE unsigned int& getMinIterations() {
			return m_minIterations;
		}
		FORCE_INLINE void setMinIterations(const unsigned int& value) {
			m_minIterations = value;
		}

		FORCE_INLINE const Real& getMaxDensityDeviation() const {
			return m_maxDensityDeviation;
		}
		FORCE_INLINE Real& getMaxDensityDeviation() {
			return m_maxDensityDeviation;
		}
		FORCE_INLINE void setMaxDensityDeviation(const Real& value) {
			m_maxDensityDeviation = value;
		}

		int getKernel() const {
			return m_kernelMethod;
		}
		void setKernel(int val);
		int getGradKernel() const {
			return m_gradKernelMethod;
		}
		void setGradKernel(int val);
	};
}