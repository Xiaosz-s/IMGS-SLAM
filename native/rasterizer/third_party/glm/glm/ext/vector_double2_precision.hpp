#pragma once
#include "../detail/type_vec2.hpp"

namespace glm
{

	typedef vec<2, double, highp>		highp_dvec2;

	typedef vec<2, double, mediump>		mediump_dvec2;

	typedef vec<2, double, lowp>		lowp_dvec2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::highp_dvec2>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_dvec2>::value);
static_assert(std::is_trivially_default_constructible<glm::lowp_dvec2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::highp_dvec2>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_dvec2>::value);
static_assert(std::is_trivially_copy_assignable<glm::lowp_dvec2>::value);
static_assert(std::is_trivially_copyable<glm::highp_dvec2>::value);
static_assert(std::is_trivially_copyable<glm::mediump_dvec2>::value);
static_assert(std::is_trivially_copyable<glm::lowp_dvec2>::value);
static_assert(std::is_copy_constructible<glm::highp_dvec2>::value);
static_assert(std::is_copy_constructible<glm::mediump_dvec2>::value);
static_assert(std::is_copy_constructible<glm::lowp_dvec2>::value);
static_assert(glm::highp_dvec2::length() == 2);
static_assert(glm::mediump_dvec2::length() == 2);
static_assert(glm::lowp_dvec2::length() == 2);
static_assert(sizeof(glm::highp_dvec2) == sizeof(glm::mediump_dvec2));
static_assert(sizeof(glm::highp_dvec2) == sizeof(glm::lowp_dvec2));

