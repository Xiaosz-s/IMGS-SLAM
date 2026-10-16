#pragma once
#include "../detail/type_vec3.hpp"

namespace glm
{

	typedef vec<3, float, defaultp>		vec3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::vec3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::vec3>::value);
static_assert(std::is_trivially_copyable<glm::vec3>::value);
static_assert(std::is_copy_constructible<glm::vec3>::value);
static_assert(glm::vec3::length() == 3);
