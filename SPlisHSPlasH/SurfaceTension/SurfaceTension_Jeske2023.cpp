#include "SurfaceTension_Jeske2023.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"

using namespace SPH;
using namespace GenParam;

int SurfaceTension_Jeske2023::ITERATIONS = -1;
int SurfaceTension_Jeske2023::MAX_ITERATIONS = -1;
int SurfaceTension_Jeske2023::MAX_ERROR = -1;
int SurfaceTension_Jeske2023::XSPH = -1;

int SurfaceTension_Jeske2023::VISCOSITY_COEFFICIENT = -1;
int SurfaceTension_Jeske2023::VISCOSITY_COEFFICIENT_BOUNDARY = -1;

SurfaceTension_Jeske2023::SurfaceTension_Jeske2023(FluidModel *model) :
        SurfaceTensionBase(model), m_vDiff(), m_gradRho() {
    m_viscosity = 0;
    m_boundaryViscosity = 0;

    m_maxIter = 100;
    m_maxError = static_cast<Real>(0.001); // This seems to be a better default than 1e-2 avoiding spinning problems etc
    m_iterations = 0;
    m_xsph = 0;

    m_vDiff.resize(model->numParticles(), Vector3r::Zero());
    m_gradRho.resize(model->numParticles(), static_cast<Real>(0.0));
    m_surfaceEnergy.resize(model->numParticles(), static_cast<Real>(0.0));
    m_color.resize(model->numParticles(), static_cast<Real>(0.0));
    m_colorGrad.resize(model->numParticles(), Vector3r::Zero());
    m_nonlinearAcc.resize(model->numParticles(), Vector3r::Zero());
    m_nonlinearRes.resize(model->numParticles(), Vector3r::Zero());
    m_nonlinearGrad.resize(model->numParticles(), Vector3r::Zero());


    model->addField( {"velocity difference", FieldType::Vector3, [&](const unsigned int i) -> Real * { return &m_vDiff[i][0]; }, true});
    model->addField( {"density gradient", FieldType::Scalar, [&](const unsigned int i) -> Real * { return &m_gradRho[i]; }, true});
    model->addField( {"surface energy", FieldType::Scalar, [&](const unsigned int i) -> Real * { return &m_surfaceEnergy[i]; }, true});
    model->addField( {"surface color", FieldType::Scalar, [&](const unsigned int i) -> Real * { return &m_color[i]; }, true});
    model->addField( {"surface color grad", FieldType::Vector3, [&](const unsigned int i) -> Real * { return &m_colorGrad[i][0]; }, true});
    model->addField( {"surface nonlinear acc", FieldType::Vector3, [&](const unsigned int i) -> Real * { return &m_nonlinearAcc[i][0]; }, true});
    model->addField( {"surface nonlinear res", FieldType::Vector3, [&](const unsigned int i) -> Real * { return &m_nonlinearRes[i][0]; }, true});
    model->addField( {"surface nonlinear grad", FieldType::Vector3, [&](const unsigned int i) -> Real * { return &m_nonlinearGrad[i][0]; }, true});
}

SurfaceTension_Jeske2023::~SurfaceTension_Jeske2023(void) {
    m_model->removeFieldByName("velocity difference");
    m_model->removeFieldByName("density gradient");
    m_model->removeFieldByName("surface energy");
    m_model->removeFieldByName("surface color");
    m_model->removeFieldByName("surface color grad");
    m_model->removeFieldByName("surface nonlinear acc");
    m_model->removeFieldByName("surface nonlinear res");
    m_model->removeFieldByName("surface nonlinear grad");

    m_vDiff.clear();
}

void SurfaceTension_Jeske2023::initParameters() {
    SurfaceTensionBase::initParameters();

    VISCOSITY_COEFFICIENT = createNumericParameter("surfaceTensionViscosity", "Viscosity coefficient ", &m_viscosity);
    setGroup(VISCOSITY_COEFFICIENT, "Fluid Model|Surface tension");
    setDescription(VISCOSITY_COEFFICIENT, "Coefficient for the viscosity force computation.");
    RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT));
    rparam->setMinValue(0.0);

    VISCOSITY_COEFFICIENT_BOUNDARY = createNumericParameter("surfaceTensionViscosityBoundary", "Viscosity coefficient (Boundary)", &m_boundaryViscosity);
    setGroup(VISCOSITY_COEFFICIENT_BOUNDARY, "Fluid Model|Surface tension");
    setDescription(VISCOSITY_COEFFICIENT_BOUNDARY, "Coefficient for the viscosity force computation at the boundary.");
    rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT_BOUNDARY));
    rparam->setMinValue(0.0);

    ITERATIONS = createNumericParameter("surfaceTensionIterations", "Iterations", &m_iterations);
    setGroup(ITERATIONS, "Fluid Model|Surface tension");
    setDescription(ITERATIONS, "Iterations required by the viscosity solver.");
    getParameter(ITERATIONS)->setReadOnly(true);

    MAX_ITERATIONS = createNumericParameter("surfaceTensionMaxIter", "Max. iterations (visco)", &m_maxIter);
    setGroup(MAX_ITERATIONS, "Fluid Model|Surface tension");
    setDescription(MAX_ITERATIONS, "Max. iterations of the viscosity solver.");
    static_cast<NumericParameter<unsigned int> *>(getParameter(MAX_ITERATIONS))->setMinValue(1);

    MAX_ERROR = createNumericParameter("surfaceTensionMaxError", "Max. surface tension error", &m_maxError);
    setGroup(MAX_ERROR, "Fluid Model|Surface tension");
    setDescription(MAX_ERROR, "Max. error of the surface tension solver.");
    // rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR));
    // rparam->setMinValue(static_cast<Real>(1e-6));

    XSPH = createNumericParameter("surfaceTensionXSPH", "XSPH Smoothing Factor", &m_xsph);
    setGroup(XSPH, "Fluid Model|Surface tension");
    setDescription(XSPH, "Factor for xsph smoothing");
}

#ifdef USE_AVX

void SurfaceTension_Jeske2023::step() {
    const int numParticles = (int)m_model->numActiveParticles();
    // prevent solver from running with a zero-length vector
    if (numParticles == 0)
        return;
    const Real density0 = m_model->getDensity0();
    const Real h = TimeManager::getCurrent()->getTimeStepSize();

    // Compute density gradient first
    computeDensityGradient();

    //////////////////////////////////////////////////////////////////////////
    // Init linear system solver and preconditioner
    //////////////////////////////////////////////////////////////////////////
    MatrixReplacement A(3 * m_model->numActiveParticles(), matrixVecProd, (void*)this);

    m_solver.setTolerance(m_maxError);
    m_solver.setMaxIterations(m_maxIter);
    m_solver.compute(A);

    VectorXr b(3 * numParticles);
    VectorXr x(3 * numParticles);
    VectorXr g(3 * numParticles);

    computeRHS(b, g);

    //////////////////////////////////////////////////////////////////////////
    // Solve linear system
    //////////////////////////////////////////////////////////////////////////
    START_TIMING("CG solve");
    x = m_solver.solveWithGuess(b, g);
    m_iterations = (int)m_solver.iterations();
    STOP_TIMING_AVG;
    INCREASE_COUNTER("Surface tension iterations", static_cast<Real>(m_iterations));

    applyForces(x);
}

void SurfaceTension_Jeske2023::matrixVecProd(const Real *vec, Real *result, void *userData) {
    Simulation *sim = Simulation::getCurrent();
    SurfaceTension_Jeske2023 *surf_tens = (SurfaceTension_Jeske2023 *) userData;
    FluidModel *model = surf_tens->getModel();
    const unsigned int numParticles = model->numActiveParticles();
    const unsigned int fluidModelIndex = model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();

    const Real h = sim->getSupportRadius();
    const Real h2 = h * h;
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();
    const Real density0 = model->getDensity0();
    const Real mu = surf_tens->m_viscosity * density0;
    const Real mub = surf_tens->m_boundaryViscosity * density0;
    const Real mu_xsph = surf_tens->m_xsph;
    const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2 * h;

    const Real cohesion = surf_tens->m_surfaceTension;
    const Real adhesion = surf_tens->m_surfaceTensionBoundary;
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    const Real W_min = sim->W(Vector3r(diameter, 0, 0));

    Real d = 10.0;
    if (sim->is2DSimulation())
        d = 8.0;

    const Scalarf8 d_mu_rho0(d * mu * density0);
    const Scalarf8 h2_001(0.01f * h2);
    const Scalarf8 W_min_avx(W_min);
    const Scalarf8 dt_avx(dt);
    const Scalarf8 diameter2_avx(diameter * diameter);
    const Scalarf8 density0_avx(density0);
    const Scalarf8 dt_cohesion_avx(dt * cohesion);
    const Scalarf8 mu_xsph_avx(mu_xsph);
    Vector3f8 zero;
    zero.setZero();

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int) numParticles; i++) {
            const Vector3r &xi = model->getPosition(i);
            Vector3r ai;
            ai.setZero();
            const Real density_i = model->getDensity(i);
            const Vector3r &vi = Eigen::Map<const Vector3r>(&vec[3 * i]);
            const Vector3r &vi_old = model->getVelocity(i);
            // const Vector3r vi_old = Vector3r::Zero();
            const Real m_i = model->getMass(i);

            const Vector3f8 xi_avx(xi);
            const Vector3f8 vi_avx(vi);
            const Vector3f8 vi_old_avx(vi_old);
            const Scalarf8 densityGrad_i_avx(surf_tens->getDensityGrad(i));
            Vector3f8 delta_ai_avx;
            const Scalarf8 density_i_avx(density_i);
            delta_ai_avx.setZero();

            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            forall_fluid_neighbors_in_same_phase_avx(
				compute_Vj(model);
				compute_Vj_gradW_samephase();
				const Scalarf8 density_j_avx = convert_one(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getDensity(0), count);
                const Scalarf8 densityGrad_j_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &surf_tens->getDensityGrad(0), count);
                const Scalarf8 Vj_avx = convert_zero(model->getVolume(0), count);
                const Scalarf8 mass_j_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getMass(0), count);
				const Vector3f8 xixj = xi_avx - xj_avx;
				const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &vec[0], count);
                const Vector3f8 vj_old_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getVelocity(0), count);

				delta_ai_avx += V_gradW * ((d_mu_rho0 / density_j_avx) * (vi_avx - vj_avx).dot(xixj) / (xixj.squaredNorm() + h2_001));
                delta_ai_avx -= (vi_avx - vj_avx) * mu_xsph_avx * mass_j_avx/density_j_avx * CubicKernel_AVX::W(xixj) * density_i_avx / dt_avx;

                const Scalarf8 W_cohesion = Scalarf8(10.f/7) * min(CubicKernel_AVX::W(xixj), W_min_avx);
                const Vector3f8 gradW = CubicKernel_AVX::gradW(xixj);
                const Vector3f8 gradW_cohesion = Vector3f8::blend((xixj).squaredNorm() > diameter2_avx, gradW, zero) * Scalarf8(10.f/7) ;
                const Scalarf8 rhoij_avx = density_i_avx + density_j_avx;

                const Scalarf8 W_firstOrder = convert_zero(2.0, count) / rhoij_avx * (W_cohesion + dt_avx * (gradW_cohesion.dot(vi_old_avx - vj_old_avx)
                                                    -  W_cohesion * (densityGrad_i_avx + densityGrad_j_avx) / rhoij_avx)
                                                    );

                delta_ai_avx -= (vi_avx - vj_avx) * dt_cohesion_avx * density_i_avx * W_firstOrder;
           
			);

            ai[0] += delta_ai_avx.x().reduce();
            ai[1] += delta_ai_avx.y().reduce();
            ai[2] += delta_ai_avx.z().reduce();
           
            //////////////////////////////////////////////////////////////////////////
            // Boundary
            //////////////////////////////////////////////////////////////////////////
            if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
            {
                forall_boundary_neighbors(
                        const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                        const Vector3r xixj = xi - xj;
                        const Vector3r gradW = sim->gradW(xixj);
                // bm_neighbor->addForce(xj, -model->getMass(i) / density_i * a);
                );
            }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
            {
                forall_density_maps(
                        const Vector3r xixj = xi - xj;
                        Vector3r normal = -xixj;
                );
            }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
            {
                // ToDo: implement adhesion for moving boundaries, i.e. vj not equal to zero
                forall_volume_maps(
                        Vector3r vj_old;
                        bm_neighbor->getPointVelocity(xj, vj_old);
                        const Vector3r xixj = xi - xj;
                        //const Real rho_iT1 = surf_tens->m_rhoT1[i];
                        const Real m_j = Vj * density0;
                        const Real density_j = m_j / Vj;
                        //const Real Vij = (m_i/rho_iT1 + Vj)*0.5;
                        const Real W_cohesion = std::min(sim->W(xixj), W_min);
                        const Vector3r gradW = sim->gradW(xixj);
                        const Vector3r gradW_cohesion = (xixj).norm() > diameter ? gradW : Vector3r::Zero();

                        // First order
                        // ai -= (vi * dt - vj * dt) * 2/(density_i + density_j) * cohesion * density_i * (std::min(sim->W(xixj), W_min) + dt * gradW_cohesion.dot(vi_old) - dt*gradW_cohesion.dot(vj_old));

                        const Real rhoij = density_i + density_j;
                        const Real m_ij = (m_i + m_j)/m_i;
                        const Real W_firstOrder = m_ij * (W_cohesion / rhoij + dt * gradW_cohesion.dot(vi_old) / rhoij
                                                            - dt * gradW_cohesion.dot(vj_old) / rhoij
                                                            - dt * W_cohesion * (surf_tens->getDensityGrad(i)
                                                            // + surf_tens->getDensityGrad(neighborIndex)
                                                            ) / (rhoij * rhoij)
                        );
                        ai -= (vi * dt
                                // - vj * dt
                                ) * adhesion * density_i * W_firstOrder;
                        //bm_neighbor->addForce(xi, -model->getMass(i) * ai);
                );
            }

            if (mub != 0.0)
            {
                if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
                {
                    forall_boundary_neighbors(
                            const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                            const Vector3r xixj = xi - xj;
                            const Vector3r gradW = sim->gradW(xixj);
                            const Vector3r a = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
                            ai += a;
                    // bm_neighbor->addForce(xj, -model->getMass(i) / density_i * a);
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
                {
                    forall_density_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = surf_tens->m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1 * dist;
                                const Vector3r x2 = xj + t1 * dist;
                                const Vector3r x3 = xj - t2 * dist;
                                const Vector3r x4 = xj + t2 * dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (vi).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (vi).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (vi).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (vi).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                ai += a1 + a2 + a3 + a4;

                                //bm_neighbor->addForce(x1, -model->getMass(i)/density_i * a1);
                                //bm_neighbor->addForce(x2, -model->getMass(i)/density_i * a2);
                                //bm_neighbor->addForce(x3, -model->getMass(i)/density_i * a3);
                                //bm_neighbor->addForce(x4, -model->getMass(i)/density_i * a4);
                            }
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
                {
                    forall_volume_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = surf_tens->m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1*dist;
                                const Vector3r x2 = xj + t1*dist;
                                const Vector3r x3 = xj - t2*dist;
                                const Vector3r x4 = xj + t2*dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * Vj;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (vi).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (vi).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (vi).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (vi).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                ai += a1 + a2 + a3 + a4;

                                //bm_neighbor->addForce(x1, -model->getMass(i)/density_i * a1);
                                //bm_neighbor->addForce(x2, -model->getMass(i)/density_i * a2);
                                //bm_neighbor->addForce(x3, -model->getMass(i)/density_i * a3);
                                //bm_neighbor->addForce(x4, -model->getMass(i)/density_i * a4);
                            }
                    );
                }
            }

            result[3 * i] = vec[3 * i] - dt / density_i * ai[0];
            result[3 * i + 1] = vec[3 * i + 1] - dt / density_i * ai[1];
            result[3 * i + 2] = vec[3 * i + 2] - dt / density_i * ai[2];
        }
    }
}

void SurfaceTension_Jeske2023::computeDensityGradient() {
    const int numParticles = (int) m_model->numActiveParticles();
    Simulation *sim = Simulation::getCurrent();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    FluidModel *model = m_model;
    const Real density0 = model->getDensity0();
    const Scalarf8 density0_avx(density0);
    Vector3f8 zero;
    zero.setZero();

    //////////////////////////////////////////////////////////////////////////
    // Compute RHS
    //////////////////////////////////////////////////////////////////////////
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static) nowait
        for (int i = 0; i < (int) numParticles; i++) {
            const Vector3r &xi = m_model->getPosition(i);
            const Vector3r &vi = m_model->getVelocity(i);
            Scalarf8 gradRho_avx(0.0);
            const Vector3f8 xi_avx(xi);
            const Vector3f8 vi_avx(vi);

            // Compute density gradient
            forall_fluid_neighbors_in_same_phase_avx(
                compute_Vj(model);
                compute_Vj_gradW_samephase();

                const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getVelocity(0), count);
                const Vector3f8 xixj = xi_avx - xj_avx;

                gradRho_avx += density0_avx * (vi_avx - vj_avx).dot(V_gradW);
                );
            m_gradRho[i] = gradRho_avx.reduce();
        }
    }
}

void SurfaceTension_Jeske2023::computeRHS(VectorXr &b, VectorXr &g) {
    const int numParticles = (int) m_model->numActiveParticles();
    Simulation *sim = Simulation::getCurrent();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    FluidModel *model = m_model;
    const Real h = sim->getSupportRadius();
    const Real h2 = h * h;
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();
    const Real density0 = m_model->getDensity0();
    const Real mu = m_viscosity * density0;
    const Real mub = m_boundaryViscosity * density0;
    const Real cohesion = m_surfaceTension;
    const Real adhesion = m_surfaceTensionBoundary;
    const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2 * h;
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    const Real W_min = sim->W(Vector3r(diameter, 0, 0));
    Real d = 10.0;
    if (sim->is2DSimulation())
        d = 8.0;

    const Scalarf8 W_min_avx(W_min);
    const Scalarf8 dt_avx(dt);
    const Scalarf8 diameter2_avx(diameter * diameter);
    const Scalarf8 density0_avx(density0);
    const Scalarf8 cohesion_avx(cohesion);
    Vector3f8 zero;
    zero.setZero();

    //////////////////////////////////////////////////////////////////////////
    // Compute RHS
    //////////////////////////////////////////////////////////////////////////
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static) nowait
        for (int i = 0; i < (int) numParticles; i++) {
            const Vector3r &vi = m_model->getVelocity(i);
            const Vector3r vi_old = vi;
            // const Vector3r vi_old = Vector3r::Zero();
            const Vector3r &xi = m_model->getPosition(i);
            const Real density_i = m_model->getDensity(i);
            const Real m_i = m_model->getMass(i);

            const Vector3f8 xi_avx(xi);
            const Vector3f8 vi_avx(vi);
            Vector3f8 delta_bi_avx;
            const Scalarf8 density_i_avx(density_i);
            delta_bi_avx.setZero();
            const Scalarf8 densityGrad_i_avx(getDensityGrad(i));

            // Implicit cohesion
            forall_fluid_neighbors_in_same_phase_avx(
                compute_Vj(model);
                compute_Vj_gradW_samephase();
                const Vector3f8 xixj = xi_avx - xj_avx;
                const Scalarf8 Vj_avx = convert_zero(model->getVolume(0), count);
                const Scalarf8 density_j_avx = convert_one(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getDensity(0), count);
                const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getVelocity(0), count);
                const Scalarf8 densityGrad_j_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &getDensityGrad(0), count);

                // Cohesion
                const Scalarf8 W_cohesion = Scalarf8(10.f/7) * min(CubicKernel_AVX::W(xixj), W_min_avx);
                const Vector3f8 gradW = CubicKernel_AVX::gradW(xixj);
                const Vector3f8 gradW_cohesion = Vector3f8::blend((xixj).squaredNorm() > diameter2_avx, gradW, zero) * Scalarf8(10.f/7);

                // First order
                const Scalarf8 rhoij_avx = density_i_avx + density_j_avx;
                const Scalarf8 W_firstOrder = convert_zero(2.0, count) / rhoij_avx * (W_cohesion + dt_avx * (gradW_cohesion.dot(vi_avx)
                                                    - gradW_cohesion.dot(vj_avx)
                                                    -  W_cohesion * (densityGrad_i_avx + densityGrad_j_avx) / rhoij_avx)
                                                    );

                delta_bi_avx += xixj * cohesion_avx * density_i_avx * W_firstOrder;
            );
            Vector3r bi = delta_bi_avx.reduce();

            //////////////////////////////////////////////////////////////////////////
            // Boundary
            //////////////////////////////////////////////////////////////////////////
            if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
            {
                forall_boundary_neighbors(
                        const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                        const Vector3r xixj = xi - xj;
                        const Vector3r gradW = sim->gradW(xixj);
                // bm_neighbor->addForce(xj, -model->getMass(i) / density_i * a);
                );
            }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
            {
                forall_density_maps(
                        const Vector3r xixj = xi - xj;
                        Vector3r normal = -xixj;
                );
            }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
            {
                // ToDo: implement adhesion for moving boundaries, i.e. vj not equal to zero
                forall_volume_maps(
                        Vector3r vj;
                        bm_neighbor->getPointVelocity(xj, vj);
                        const Vector3r xixj = xi - xj;
                        //const Real rho_iT1 = m_rhoT1[i];
                        const Real m_j = Vj * density0;
                        const Real density_j = m_j / Vj;
                        //const Real Vij = (m_i/rho_iT1 + Vj)*0.5;
                        const Real W_cohesion = std::min(sim->W(xixj), W_min);
                        const Vector3r gradW = sim->gradW(xixj);
                        const Vector3r gradW_cohesion = (xixj).norm() > diameter ? gradW : Vector3r::Zero();

                        const Real m_ij = (m_i + m_j)/m_i;
                        const Real rhoij = density_i + density_j;
                        // TODO: some issue with adhesion here
                        const Real W_firstOrder = m_ij * ( W_cohesion / rhoij
                                                            + dt * gradW_cohesion.dot(vi) / rhoij
                                                            - dt * gradW_cohesion.dot(vj) / rhoij
                                                            - dt * W_cohesion * (getDensityGrad(i)
                                                            // + getDensityGrad(neighborIndex)
                                                            ) / (rhoij * rhoij)
                        );
                        const Vector3r a = adhesion * xixj * density_i * W_firstOrder;
                        bi += a;
                        bm_neighbor->addForce(xi, -m_model->getMass(i) / density_i * a);
                );
            }
            if (mub != 0.0)
            {
                if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
                {
                    forall_boundary_neighbors(
                            const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                            const Vector3r xixj = xi - xj;
                            const Vector3r gradW = sim->gradW(xixj);
                            const Vector3r a = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
                            bi += a;
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
                {
                    forall_density_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1 * dist;
                                const Vector3r x2 = xj + t1 * dist;
                                const Vector3r x3 = xj - t2 * dist;
                                const Vector3r x4 = xj + t2 * dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                bi += a1 + a2 + a3 + a4;
                            }
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
                {
                    forall_volume_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1*dist;
                                const Vector3r x2 = xj + t1*dist;
                                const Vector3r x3 = xj - t2*dist;
                                const Vector3r x4 = xj + t2*dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * Vj;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                bi += a1 + a2 + a3 + a4;
                            }
                    );
                }
            }

            b[3 * i] = vi[0] - dt / density_i * bi[0];
            b[3 * i + 1] = vi[1] - dt / density_i * bi[1];
            b[3 * i + 2] = vi[2] - dt / density_i * bi[2];

            // Warmstart
            g[3 * i] = vi[0] + m_vDiff[i][0];
            g[3 * i + 1] = vi[1] + m_vDiff[i][1];
            g[3 * i + 2] = vi[2] + m_vDiff[i][2];
        }
    }
}


#else 

void SurfaceTension_Jeske2023::matrixVecProd(const Real *vec, Real *result, void *userData) {
    Simulation *sim = Simulation::getCurrent();
    SurfaceTension_Jeske2023 *surf_tens = (SurfaceTension_Jeske2023 *) userData;
    FluidModel *model = surf_tens->getModel();
    const unsigned int numParticles = model->numActiveParticles();
    const unsigned int fluidModelIndex = model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();

    const Real h = sim->getSupportRadius();
    const Real h2 = h * h;
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();
    const Real density0 = model->getDensity0();
    const Real mu = surf_tens->m_viscosity * density0;
    const Real mub = surf_tens->m_boundaryViscosity * density0;
    const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2 * h;

    const Real cohesion = surf_tens->m_surfaceTension;
    const Real adhesion = surf_tens->m_surfaceTensionBoundary;
    const Real mu_xsph = surf_tens->m_xsph;
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    const Real W_min = sim->W(Vector3r(diameter, 0, 0));

    Real d = 10.0;
    if (sim->is2DSimulation())
        d = 8.0;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int) numParticles; i++) {
            const Vector3r &xi = model->getPosition(i);
            Vector3r ai;
            ai.setZero();
            const Real density_i = model->getDensity(i);
            const Vector3r &vi = Eigen::Map<const Vector3r>(&vec[3 * i]);
            const Vector3r &vi_old = model->getVelocity(i);
            // const Vector3r vi_old = Vector3r::Zero();
            const Real m_i = model->getMass(i);

            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            forall_fluid_neighbors_in_same_phase(
                    const Real density_j = model->getDensity(neighborIndex);
                    const Real mass_j = model->getMass(neighborIndex);
                    const Vector3r gradW = sim->gradW(xi - xj);

                    const Vector3r &vj = Eigen::Map<const Vector3r>(&vec[3 * neighborIndex]);
                    const Vector3r &vj_old = model->getVelocity(neighborIndex);
                    // const Vector3r vj_old = Vector3r::Zero();
                    const Vector3r xixj = xi - xj;
                    const Real m_ij = m_i + mass_j;

                    ai += d * mu * (model->getMass(neighborIndex) / density_j) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
                    ai -= (vi - vj) * mu_xsph * mass_j / density_j * CubicKernel::W(xixj) * density_i / dt;

                    // Cohesion
                    const Real W_cohesion = (10./7.) * std::min(sim->W(xixj), W_min);
                    Vector3r gradW_cohesion = (xi - xj).norm() > diameter ? gradW : Vector3r::Zero();
                    gradW_cohesion *= (10./7.);

                    // First order
                    // ai -= (vi * dt - vj * dt) * 2/(density_i + density_j) * cohesion * density_i * (std::min(sim->W(xixj), W_min) + dt * gradW_cohesion.dot(vi_old) - dt*gradW_cohesion.dot(vj_old));

                    const Real rhoij = density_i + density_j;
                    const Real W_firstOrder = 2.0 * (W_cohesion / rhoij + dt * gradW_cohesion.dot(vi_old) / rhoij // - dt * W_cohesion * (surf_tens->getDensityGrad(i) - m_i * gradW ).dot(vi_old) / (rhoij * rhoij)
                                                    - dt * gradW_cohesion.dot(vj_old) / rhoij // - dt * W_cohesion * (surf_tens->getDensityGrad(neighborIndex) - mass_j * gradW ).dot(vj_old) / (rhoij * rhoij)
                                                    - dt * W_cohesion * (surf_tens->getDensityGrad(i) + surf_tens->getDensityGrad(neighborIndex)) / (rhoij * rhoij)
                                                    );
                    ai -= (vi * dt - vj * dt) * cohesion * density_i * W_firstOrder;

            );

            //////////////////////////////////////////////////////////////////////////
            // Boundary
            //////////////////////////////////////////////////////////////////////////
            if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
            {
                forall_boundary_neighbors(
                        const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                        const Vector3r xixj = xi - xj;
                        const Vector3r gradW = sim->gradW(xixj);
                // bm_neighbor->addForce(xj, -model->getMass(i) / density_i * a);
                );
            }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
            {
                forall_density_maps(
                        const Vector3r xixj = xi - xj;
                        Vector3r normal = -xixj;
                );
            }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
            {
                // ToDo: implement adhesion for moving boundaries, i.e. vj not equal to zero
                forall_volume_maps(
                        Vector3r vj_old;
                        bm_neighbor->getPointVelocity(xj, vj_old);
                        const Vector3r xixj = xi - xj;
                        const Real m_j = Vj * density0;
                        const Real density_j = m_j / Vj;
                        const Real W_cohesion = std::min(sim->W(xixj), W_min);
                        const Vector3r gradW = sim->gradW(xixj);
                        const Vector3r gradW_cohesion = (xixj).norm() > diameter ? gradW : Vector3r::Zero();

                        const Real rhoij = density_i + density_j;
                        const Real m_ij = (m_i + m_j)/m_i;
                        const Real W_firstOrder = m_ij * (W_cohesion / rhoij + dt * gradW_cohesion.dot(vi_old) / rhoij
                                                            - dt * gradW_cohesion.dot(vj_old) / rhoij
                                                            - dt * W_cohesion * (surf_tens->getDensityGrad(i)
                                                            // + surf_tens->getDensityGrad(neighborIndex)
                                                            ) / (rhoij * rhoij)
                        );
                        ai -= (vi * dt
                                ) * adhesion * density_i * W_firstOrder;
                );
            }

            if (mub != 0.0)
            {
                if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
                {
                    forall_boundary_neighbors(
                            const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                            const Vector3r xixj = xi - xj;
                            const Vector3r gradW = sim->gradW(xixj);
                            const Vector3r a = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
                            ai += a;
                    // bm_neighbor->addForce(xj, -model->getMass(i) / density_i * a);
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
                {
                    forall_density_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = surf_tens->m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1 * dist;
                                const Vector3r x2 = xj + t1 * dist;
                                const Vector3r x3 = xj - t2 * dist;
                                const Vector3r x4 = xj + t2 * dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (vi).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (vi).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (vi).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (vi).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                ai += a1 + a2 + a3 + a4;

                                // bm_neighbor->addForce(x1, -model->getMass(i)/density_i * a1);
                                // bm_neighbor->addForce(x2, -model->getMass(i)/density_i * a2);
                                // bm_neighbor->addForce(x3, -model->getMass(i)/density_i * a3);
                                // bm_neighbor->addForce(x4, -model->getMass(i)/density_i * a4);
                            }
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
                {
                    forall_volume_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = surf_tens->m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1*dist;
                                const Vector3r x2 = xj + t1*dist;
                                const Vector3r x3 = xj - t2*dist;
                                const Vector3r x4 = xj + t2*dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * Vj;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (vi).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (vi).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (vi).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (vi).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                ai += a1 + a2 + a3 + a4;

                                // bm_neighbor->addForce(x1, -model->getMass(i)/density_i * a1);
                                // bm_neighbor->addForce(x2, -model->getMass(i)/density_i * a2);
                                // bm_neighbor->addForce(x3, -model->getMass(i)/density_i * a3);
                                // bm_neighbor->addForce(x4, -model->getMass(i)/density_i * a4);
                            }
                    );
                }
            }

            result[3 * i] = vec[3 * i] - dt / density_i * ai[0];
            result[3 * i + 1] = vec[3 * i + 1] - dt / density_i * ai[1];
            result[3 * i + 2] = vec[3 * i + 2] - dt / density_i * ai[2];
        }
    }
}

void SurfaceTension_Jeske2023::step() {
    const int numParticles = (int)m_model->numActiveParticles();
    // prevent solver from running with a zero-length vector
    if (numParticles == 0)
        return;
    const Real density0 = m_model->getDensity0();
    const Real h = TimeManager::getCurrent()->getTimeStepSize();

    // Compute xsph first
    // step_xsph();
//    solve_nonlinear();

    // Compute density gradient first
    computeDensityGradient();

    //////////////////////////////////////////////////////////////////////////
    // Init linear system solver and preconditioner
    //////////////////////////////////////////////////////////////////////////
    MatrixReplacement A(3 * m_model->numActiveParticles(), matrixVecProd, (void*)this);
    // m_solver.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElement, (void*)this);

    m_solver.setTolerance(m_maxError);
    m_solver.setMaxIterations(m_maxIter);
    m_solver.compute(A);

    VectorXr b(3 * numParticles);
    VectorXr x(3 * numParticles);
    VectorXr g(3 * numParticles);

    computeRHS(b, g);

    //////////////////////////////////////////////////////////////////////////
    // Solve linear system
    //////////////////////////////////////////////////////////////////////////
    START_TIMING("CG solve");
    x = m_solver.solveWithGuess(b, g);
    m_iterations = (int)m_solver.iterations();
    STOP_TIMING_AVG;
    INCREASE_COUNTER("Surface tension iterations", static_cast<Real>(m_iterations));

    applyForces(x);

}

void SurfaceTension_Jeske2023::computeDensityGradient() {
    const int numParticles = (int) m_model->numActiveParticles();
    Simulation *sim = Simulation::getCurrent();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    FluidModel *model = m_model;
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    // Reset density gradients
    std::fill(m_gradRho.begin(), m_gradRho.end(), static_cast<Real>(0.0));

    //////////////////////////////////////////////////////////////////////////
    // Compute RHS
    //////////////////////////////////////////////////////////////////////////
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static) nowait
        for (int i = 0; i < (int) numParticles; i++) {
            const Vector3r &xi = m_model->getPosition(i);
            const Vector3r &vi = m_model->getVelocity(i);

            // Compute density gradient
            forall_fluid_neighbors_in_same_phase(
                    const Vector3r xixj = xi - xj;
                    const Vector3r &vj = m_model->getVelocity(neighborIndex);
                    const Vector3r gradW = sim->gradW(xixj);
                    const Real m_j = model->getMass(neighborIndex);
                    m_gradRho[i] += m_j * (vi - vj).dot(gradW);
            )
        }
    }
}

void SurfaceTension_Jeske2023::computeRHS(VectorXr &b, VectorXr &g) {
    const int numParticles = (int) m_model->numActiveParticles();
    Simulation *sim = Simulation::getCurrent();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    FluidModel *model = m_model;
    const Real h = sim->getSupportRadius();
    const Real h2 = h * h;
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();
    const Real density0 = m_model->getDensity0();
    const Real mu = m_viscosity * density0;
    const Real mub = m_boundaryViscosity * density0;
    const Real cohesion = m_surfaceTension;
    const Real adhesion = m_surfaceTensionBoundary;
    const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2 * h;
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    const Real W_min = sim->W(Vector3r(diameter, 0, 0));
    Real d = 10.0;
    if (sim->is2DSimulation())
        d = 8.0;

    //////////////////////////////////////////////////////////////////////////
    // Compute RHS
    //////////////////////////////////////////////////////////////////////////
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static) nowait
        for (int i = 0; i < (int) numParticles; i++) {
            const Vector3r &vi = m_model->getVelocity(i);
            const Vector3r vi_old = vi;
            const Vector3r &xi = m_model->getPosition(i);
            const Real density_i = m_model->getDensity(i);
            const Real m_i = m_model->getMass(i);
            Vector3r bi = Vector3r::Zero();

            // Implicit cohesion
            forall_fluid_neighbors_in_same_phase(
                    const Vector3r xixj = xi - xj;
                    const Real m_j = model->getMass(neighborIndex);
                    const Real m_ij = m_i + m_j;
                    const Real density_j = m_model->getDensity(neighborIndex);
                    const Vector3r &vj = model->getVelocity(neighborIndex);
                    const Vector3r vj_old = vj;
                    const Vector3r gradW = sim->gradW(xi - xj);

                    // Cohesion
                    const Real W_cohesion = (10./7.) * std::min(sim->W(xixj), W_min);
                    Vector3r gradW_cohesion = (xi - xj).norm() > diameter ? gradW : Vector3r::Zero();
                    gradW_cohesion *= (10./7.);

                    const Real rhoij = density_i + density_j;
                    const Real W_firstOrder = 2.0 * ( W_cohesion / rhoij
                                                    + dt * gradW_cohesion.dot(vi) / rhoij 
                                                    - dt * gradW_cohesion.dot(vj) / rhoij 
                                                    - dt * W_cohesion * (getDensityGrad(i) + getDensityGrad(neighborIndex)) / (rhoij * rhoij)
                    );
                    bi += cohesion * xixj * density_i * W_firstOrder;
                        
            )

            //////////////////////////////////////////////////////////////////////////
            // Boundary
            //////////////////////////////////////////////////////////////////////////
            if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
            {
                forall_boundary_neighbors(
                        const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                        const Vector3r xixj = xi - xj;
                        const Vector3r gradW = sim->gradW(xixj);
                // bm_neighbor->addForce(xj, -model->getMass(i) / density_i * a);
                );
            }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
            {
                forall_density_maps(
                        const Vector3r xixj = xi - xj;
                        Vector3r normal = -xixj;
                );
            }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
            {
                // ToDo: implement adhesion for moving boundaries, i.e. vj not equal to zero
                forall_volume_maps(
                        Vector3r vj;
                        bm_neighbor->getPointVelocity(xj, vj);
                        const Vector3r xixj = xi - xj;
                        //const Real rho_iT1 = m_rhoT1[i];
                        const Real m_j = Vj * density0;
                        const Real density_j = m_j / Vj;
                        //const Real Vij = (m_i/rho_iT1 + Vj)*0.5;
                        const Real W_cohesion = std::min(sim->W(xixj), W_min);
                        const Vector3r gradW = sim->gradW(xixj);
                        const Vector3r gradW_cohesion = (xixj).norm() > diameter ? gradW : Vector3r::Zero();

                        const Real m_ij = (m_i + m_j)/m_i;
                        const Real rhoij = density_i + density_j;
                        // TODO: some issue with adhesion here
                        const Real W_firstOrder = m_ij * ( W_cohesion / rhoij
                                                            + dt * gradW_cohesion.dot(vi) / rhoij
                                                            - dt * gradW_cohesion.dot(vj) / rhoij
                                                            - dt * W_cohesion * (getDensityGrad(i)
                                                            // + getDensityGrad(neighborIndex)
                                                            ) / (rhoij * rhoij)
                        );
                        bi += adhesion * xixj * density_i * W_firstOrder;

                );
            }
            if (mub != 0.0)
            {
                if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
                {
                    forall_boundary_neighbors(
                            const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                            const Vector3r xixj = xi - xj;
                            const Vector3r gradW = sim->gradW(xixj);
                            const Vector3r a = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
                            bi += a;
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
                {
                    forall_density_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1 * dist;
                                const Vector3r x2 = xj + t1 * dist;
                                const Vector3r x3 = xj - t2 * dist;
                                const Vector3r x4 = xj + t2 * dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                bi += a1 + a2 + a3 + a4;
                            }
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
                {
                    forall_volume_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1*dist;
                                const Vector3r x2 = xj + t1*dist;
                                const Vector3r x3 = xj - t2*dist;
                                const Vector3r x4 = xj + t2*dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * Vj;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                bi += a1 + a2 + a3 + a4;
                            }
                    );
                }
            }

            b[3 * i] = vi[0] - dt / density_i * bi[0];
            b[3 * i + 1] = vi[1] - dt / density_i * bi[1];
            b[3 * i + 2] = vi[2] - dt / density_i * bi[2];

            // Warmstart
            g[3 * i] = vi[0] + m_vDiff[i][0];
            g[3 * i + 1] = vi[1] + m_vDiff[i][1];
            g[3 * i + 2] = vi[2] + m_vDiff[i][2];
        }
    }
}


#endif


void SurfaceTension_Jeske2023::applyForces(const VectorXr &x) {
    const int numParticles = (int) m_model->numActiveParticles();
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int) numParticles; i++) {
            // Compute the acceleration from the velocity change
            Vector3r &ai = m_model->getAcceleration(i);
            const Vector3r newVi(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
            ai += (1.0 / dt) * (newVi - m_model->getVelocity(i));
            m_vDiff[i] = (newVi - m_model->getVelocity(i));
        }
    }
}



void SurfaceTension_Jeske2023::reset() {
    std::fill(m_vDiff.begin(), m_vDiff.end(), Vector3r::Zero());
}

void SurfaceTension_Jeske2023::performNeighborhoodSearchSort() {
     const unsigned int numPart = m_model->numActiveParticles();
     if (numPart == 0)
         return;

     Simulation *sim = Simulation::getCurrent();
     auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
     d.sort_field(&m_vDiff[0]);
}
