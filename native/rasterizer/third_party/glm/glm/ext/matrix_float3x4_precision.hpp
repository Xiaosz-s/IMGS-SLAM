#pragma once
#include "../detail/type_mat3x4.hpp"

namespace glm
{

	typedef mat<3, 4, float, lowp>		lowp_mat3x4;

	typedef mat<3, 4, float, mediump>	mediump_mat3x4;

	typedef mat<3, 4, float, highp>		highp_mat3x4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::lowp_mat3x4>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_mat3x4>::value);
static_assert(std::is_trivially_default_constructible<glm::highp_mat3x4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::lowp_mat3x4>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_mat3x4>::value);
static_assert(std::is_trivially_copy_assignable<glm::highp_mat3x4>::value);
static_assert(std::is_trivially_copyable<glm::lowp_mat3x4>::value);
static_assert(std::is_trivially_copyable<glm::mediump_mat3x4>::value);
static_assert(std::is_trivially_copyable<glm::highp_mat3x4>::value);
static_assert(std::is_copy_constructible<glm::lowp_mat3x4>::value);
static_assert(std::is_copy_constructible<glm::mediump_mat3x4>::value);
static_assert(std::is_copy_constructible<glm::highp_mat3x4>::value);
static_assert(glm::lowp_mat3x4::length() == 3);
static_assert(glm::mediump_mat3x4::length() == 3);
static_assert(glm::highp_mat3x4::length() == 3);
