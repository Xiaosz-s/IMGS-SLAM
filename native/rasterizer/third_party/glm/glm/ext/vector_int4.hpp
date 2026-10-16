#pragma once
#include "../detail/type_vec4.hpp"

namespace glm
{

	typedef vec<4, int, defaultp>		ivec4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::ivec4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::ivec4>::value);
static_assert(std::is_trivially_copyable<glm::ivec4>::value);
static_assert(std::is_copy_constructible<glm::ivec4>::value);
static_assert(glm::ivec4::length() == 4);
