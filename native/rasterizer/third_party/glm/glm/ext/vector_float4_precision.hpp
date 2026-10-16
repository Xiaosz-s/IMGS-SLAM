#pragma once
#include "../detail/type_vec4.hpp"

namespace glm
{

	typedef vec<4, float, highp>		highp_vec4;

	typedef vec<4, float, mediump>		mediump_vec4;

	typedef vec<4, float, lowp>			lowp_vec4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::highp_vec4>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_vec4>::value);
static_assert(std::is_trivially_default_constructible<glm::lowp_vec4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::highp_vec4>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_vec4>::value);
static_assert(std::is_trivially_copy_assignable<glm::lowp_vec4>::value);
static_assert(std::is_trivially_copyable<glm::highp_vec4>::value);
static_assert(std::is_trivially_copyable<glm::mediump_vec4>::value);
static_assert(std::is_trivially_copyable<glm::lowp_vec4>::value);
static_assert(std::is_copy_constructible<glm::highp_vec4>::value);
static_assert(std::is_copy_constructible<glm::mediump_vec4>::value);
static_assert(std::is_copy_constructible<glm::lowp_vec4>::value);
static_assert(glm::highp_vec4::length() == 4);
static_assert(glm::mediump_vec4::length() == 4);
static_assert(glm::lowp_vec4::length() == 4);
static_assert(sizeof(glm::highp_vec4) == sizeof(glm::mediump_vec4));
static_assert(sizeof(glm::highp_vec4) == sizeof(glm::lowp_vec4));
