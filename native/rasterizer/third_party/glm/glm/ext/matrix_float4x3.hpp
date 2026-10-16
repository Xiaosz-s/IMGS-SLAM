#pragma once
#include "../detail/type_mat4x3.hpp"

namespace glm
{

	typedef mat<4, 3, float, defaultp>			mat4x3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::mat4x3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::mat4x3>::value);
static_assert(std::is_trivially_copyable<glm::mat4x3>::value);
static_assert(std::is_copy_constructible<glm::mat4x3>::value);
static_assert(glm::mat4x3::length() == 4);
