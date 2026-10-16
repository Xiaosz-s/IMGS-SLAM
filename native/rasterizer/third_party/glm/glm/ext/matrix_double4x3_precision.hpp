#pragma once
#include "../detail/type_mat4x3.hpp"

namespace glm
{

	typedef mat<4, 3, double, lowp>		lowp_dmat4x3;

	typedef mat<4, 3, double, mediump>	mediump_dmat4x3;

	typedef mat<4, 3, double, highp>	highp_dmat4x3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::lowp_dmat4x3>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_dmat4x3>::value);
static_assert(std::is_trivially_default_constructible<glm::highp_dmat4x3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::lowp_dmat4x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_dmat4x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::highp_dmat4x3>::value);
static_assert(std::is_trivially_copyable<glm::lowp_dmat4x3>::value);
static_assert(std::is_trivially_copyable<glm::mediump_dmat4x3>::value);
static_assert(std::is_trivially_copyable<glm::highp_dmat4x3>::value);
static_assert(std::is_copy_constructible<glm::lowp_dmat4x3>::value);
static_assert(std::is_copy_constructible<glm::mediump_dmat4x3>::value);
static_assert(std::is_copy_constructible<glm::highp_dmat4x3>::value);
static_assert(glm::lowp_dmat4x3::length() == 4);
static_assert(glm::mediump_dmat4x3::length() == 4);
static_assert(glm::highp_dmat4x3::length() == 4);
