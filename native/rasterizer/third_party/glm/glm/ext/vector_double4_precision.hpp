#pragma once
#include "../detail/setup.hpp"
#include "../detail/type_vec4.hpp"

namespace glm
{

	typedef vec<4, double, highp>		highp_dvec4;

	typedef vec<4, double, mediump>		mediump_dvec4;

	typedef vec<4, double, lowp>		lowp_dvec4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::highp_dvec4>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_dvec4>::value);
static_assert(std::is_trivially_default_constructible<glm::lowp_dvec4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::highp_dvec4>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_dvec4>::value);
static_assert(std::is_trivially_copy_assignable<glm::lowp_dvec4>::value);
static_assert(std::is_trivially_copyable<glm::highp_dvec4>::value);
static_assert(std::is_trivially_copyable<glm::mediump_dvec4>::value);
static_assert(std::is_trivially_copyable<glm::lowp_dvec4>::value);
static_assert(std::is_copy_constructible<glm::highp_dvec4>::value);
static_assert(std::is_copy_constructible<glm::mediump_dvec4>::value);
static_assert(std::is_copy_constructible<glm::lowp_dvec4>::value);
static_assert(glm::highp_dvec4::length() == 4);
static_assert(glm::mediump_dvec4::length() == 4);
static_assert(glm::lowp_dvec4::length() == 4);
static_assert(sizeof(glm::highp_dvec4) == sizeof(glm::mediump_dvec4));
static_assert(sizeof(glm::highp_dvec4) == sizeof(glm::lowp_dvec4));
