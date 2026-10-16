#pragma once
#include "../detail/type_mat3x4.hpp"

namespace glm
{

	typedef mat<3, 4, float, defaultp>			mat3x4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::mat3x4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::mat3x4>::value);
static_assert(std::is_trivially_copyable<glm::mat3x4>::value);
static_assert(std::is_copy_constructible<glm::mat3x4>::value);
static_assert(glm::mat3x4::length() == 3);
