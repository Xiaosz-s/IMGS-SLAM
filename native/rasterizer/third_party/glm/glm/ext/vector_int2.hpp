#pragma once
#include "../detail/type_vec2.hpp"

namespace glm
{

	typedef vec<2, int, defaultp>		ivec2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::ivec2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::ivec2>::value);
static_assert(std::is_trivially_copyable<glm::ivec2>::value);
static_assert(std::is_copy_constructible<glm::ivec2>::value);
static_assert(glm::ivec2::length() == 2);

