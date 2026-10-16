#pragma once
#include "../detail/type_vec4.hpp"

namespace glm
{

	typedef vec<4, unsigned int, defaultp>		uvec4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::uvec4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::uvec4>::value);
static_assert(std::is_trivially_copyable<glm::uvec4>::value);
static_assert(std::is_copy_constructible<glm::uvec4>::value);
static_assert(glm::uvec4::length() == 4);

