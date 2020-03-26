#include "common.h"
#include <SPlisHSPlasH/Common.h>
#include <SPlisHSPlasH/SPHKernels.h>
#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

template<class T>
py::class_<T> define_kernel(py::module m_sub, const char* name) {
	auto cl = py::class_<T>(m_sub, name)
		.def(py::init<>())
		.def_static("getRadius", &T::getRadius)
		.def_static("setRadius", &T::setRadius)
		.def_static("W", overload_cast_<const Real>()(&T::W))
		.def_static("W", overload_cast_<const Vector3r&>()(&T::W))
		.def_static("gradW", &T::gradW)
		.def_static("W_zero", &T::W_zero);

	return cl;
}

template<class T>
py::class_<T> define_kernel_no_grad(py::module m_sub, const char* name) {
	auto cl = py::class_<T>(m_sub, name)
		.def(py::init<>())
		.def_static("getRadius", &T::getRadius)
		.def_static("setRadius", &T::setRadius)
		.def_static("W", overload_cast_<const Real>()(&T::W))
		.def_static("W", overload_cast_<const Vector3r&>()(&T::W))
		.def_static("W_zero", &T::W_zero);

	return cl;
}

void SPHKernelsModule(py::module m)
{
	auto m_sub = m.def_submodule("SPHKernels");
	auto m_sub_sub = m_sub.def_submodule("Precomputed");

	// Cubic Kernel
	auto kernel = define_kernel<SPH::CubicKernel>(m_sub, "CubicKernel");
	kernel = define_kernel<SPH::PrecomputedKernel<SPH::CubicKernel>>(m_sub_sub, "CubicKernel");

	// Cubic Kernel 2D
	kernel = define_kernel<SPH::CubicKernel2D>(m_sub, "CubicKernel2D");
	kernel = define_kernel<SPH::PrecomputedKernel<SPH::CubicKernel2D>>(m_sub_sub, "CubicKernel2D");

	// Poly6 Kernel
	kernel = define_kernel<SPH::Poly6Kernel>(m_sub, "Poly6Kernel");
	kernel.def_static("laplacianW", &SPH::Poly6Kernel::laplacianW);
	kernel = define_kernel<SPH::PrecomputedKernel<SPH::Poly6Kernel>>(m_sub_sub, "Poly6Kernel");

	// Spiky Kernel	
	kernel = define_kernel<SPH::SpikyKernel>(m_sub, "SpikyKernel");
	kernel = define_kernel<SPH::PrecomputedKernel<SPH::SpikyKernel>>(m_sub_sub, "SpikyKernel");

	// WendlandQuinticC2 Kernel	
	kernel = define_kernel<SPH::WendlandQuinticC2Kernel>(m_sub, "WendlandQuinticC2Kernel");
	kernel = define_kernel<SPH::PrecomputedKernel<SPH::WendlandQuinticC2Kernel>>(m_sub_sub, "WendlandQuinticC2Kernel");

	// WendlandQuinticC2 Kernel	2D
	kernel = define_kernel<SPH::WendlandQuinticC2Kernel2D>(m_sub, "WendlandQuinticC2Kernel2D");
	kernel = define_kernel<SPH::PrecomputedKernel<SPH::WendlandQuinticC2Kernel2D>>(m_sub_sub, "WendlandQuinticC2Kernel2D");

	// Cohesion Kernel	
	kernel = define_kernel_no_grad<SPH::CohesionKernel>(m_sub, "CohesionKernel");

	// Adhesion Kernel	
	kernel = define_kernel_no_grad<SPH::AdhesionKernel>(m_sub, "AdhesionKernel");
}