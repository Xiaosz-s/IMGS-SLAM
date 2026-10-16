#pragma once
#include "../detail/type_vec2.hpp"

namespace glm
{

	typedef vec<2, unsigned int, defaultp>		uvec2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::uvec2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::uvec2>::value);
static_assert(std::is_trivially_copyable<glm::uvec2>::value);
static_assert(std::is_copy_constructible<glm::uvec2>::value);
static_assert(glm::uvec2::length() == 2);
