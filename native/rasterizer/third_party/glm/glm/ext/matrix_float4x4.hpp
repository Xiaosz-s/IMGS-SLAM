#pragma once
#include "../detail/type_mat4x4.hpp"

namespace glm
{

	typedef mat<4, 4, float, defaultp>			mat4x4;

	typedef mat<4, 4, float, defaultp>			mat4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::mat4x4>::value);
static_assert(std::is_trivially_default_constructible<glm::mat4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::mat4x4>::value);
static_assert(std::is_trivially_copy_assignable<glm::mat4>::value);
static_assert(std::is_trivially_copyable<glm::mat4x4>::value);
static_assert(std::is_trivially_copyable<glm::mat4>::value);
static_assert(std::is_copy_constructible<glm::mat4x4>::value);
static_assert(std::is_copy_constructible<glm::mat4>::value);
static_assert(glm::mat4x4::length() == 4);
static_assert(glm::mat4::length() == 4);
