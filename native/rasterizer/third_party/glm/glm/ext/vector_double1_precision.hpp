#pragma once

#include "../detail/type_vec1.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_vector_double1_precision extension included")
#endif

namespace glm
{

	typedef vec<1, double, highp>		highp_dvec1;

	typedef vec<1, double, mediump>		mediump_dvec1;

	typedef vec<1, double, lowp>		lowp_dvec1;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::highp_dvec1>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_dvec1>::value);
static_assert(std::is_trivially_default_constructible<glm::lowp_dvec1>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::highp_dvec1>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_dvec1>::value);
static_assert(std::is_trivially_copy_assignable<glm::lowp_dvec1>::value);
static_assert(std::is_trivially_copyable<glm::highp_dvec1>::value);
static_assert(std::is_trivially_copyable<glm::mediump_dvec1>::value);
static_assert(std::is_trivially_copyable<glm::lowp_dvec1>::value);
static_assert(std::is_copy_constructible<glm::highp_dvec1>::value);
static_assert(std::is_copy_constructible<glm::mediump_dvec1>::value);
static_assert(std::is_copy_constructible<glm::lowp_dvec1>::value);
static_assert(glm::highp_dvec1::length() == 1);
static_assert(glm::mediump_dvec1::length() == 1);
static_assert(glm::lowp_dvec1::length() == 1);
static_assert(sizeof(glm::highp_dvec1) == sizeof(glm::mediump_dvec1));
static_assert(sizeof(glm::highp_dvec1) == sizeof(glm::lowp_dvec1));
