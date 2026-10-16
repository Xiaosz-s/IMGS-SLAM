#pragma once
#include "../detail/type_vec3.hpp"

namespace glm
{

	typedef vec<3, int, defaultp>		ivec3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::ivec3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::ivec3>::value);
static_assert(std::is_trivially_copyable<glm::ivec3>::value);
static_assert(std::is_copy_constructible<glm::ivec3>::value);
static_assert(glm::ivec3::length() == 3);

