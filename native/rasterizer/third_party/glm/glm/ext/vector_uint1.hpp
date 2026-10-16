#pragma once

#include "../detail/type_vec1.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_vector_uint1 extension included")
#endif

namespace glm
{

	typedef vec<1, unsigned int, defaultp>			uvec1;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::uvec1>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::uvec1>::value);
static_assert(std::is_trivially_copyable<glm::uvec1>::value);
static_assert(std::is_copy_constructible<glm::uvec1>::value);
static_assert(glm::uvec1::length() == 1);
