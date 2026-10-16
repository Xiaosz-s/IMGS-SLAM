#pragma once
#include "../detail/type_mat3x3.hpp"

namespace glm
{

	typedef mat<3, 3, float, defaultp>			mat3x3;

	typedef mat<3, 3, float, defaultp>			mat3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::mat3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::mat3>::value);
static_assert(std::is_trivially_copyable<glm::mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::mat3>::value);
static_assert(std::is_copy_constructible<glm::mat3x3>::value);
static_assert(std::is_copy_constructible<glm::mat3>::value);
static_assert(glm::mat3x3::length() == 3);
static_assert(glm::mat3::length() == 3);
