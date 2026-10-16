#pragma once
#include "../detail/type_mat2x2.hpp"

namespace glm
{

	typedef mat<4, 2, float, lowp>		lowp_mat4x2;

	typedef mat<4, 2, float, mediump>	mediump_mat4x2;

	typedef mat<4, 2, float, highp>		highp_mat4x2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::lowp_mat4x2>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_mat4x2>::value);
static_assert(std::is_trivially_default_constructible<glm::highp_mat4x2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::lowp_mat4x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_mat4x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::highp_mat4x2>::value);
static_assert(std::is_trivially_copyable<glm::lowp_mat4x2>::value);
static_assert(std::is_trivially_copyable<glm::mediump_mat4x2>::value);
static_assert(std::is_trivially_copyable<glm::highp_mat4x2>::value);
static_assert(std::is_copy_constructible<glm::lowp_mat4x2>::value);
static_assert(std::is_copy_constructible<glm::mediump_mat4x2>::value);
static_assert(std::is_copy_constructible<glm::highp_mat4x2>::value);
static_assert(glm::lowp_mat4x2::length() == 4);
static_assert(glm::mediump_mat4x2::length() == 4);
static_assert(glm::highp_mat4x2::length() == 4);
