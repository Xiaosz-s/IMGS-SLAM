#pragma once

#include "../ext/vector_int2.hpp"
#include "../ext/scalar_int_sized.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_vector_int2_sized extension included")
#endif

namespace glm
{

	typedef vec<2, int8, defaultp>		i8vec2;

	typedef vec<2, int16, defaultp>		i16vec2;

	typedef vec<2, int32, defaultp>		i32vec2;

	typedef vec<2, int64, defaultp>		i64vec2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::i8vec2>::value);
static_assert(std::is_trivially_default_constructible<glm::i16vec2>::value);
static_assert(std::is_trivially_default_constructible<glm::i32vec2>::value);
static_assert(std::is_trivially_default_constructible<glm::i64vec2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::i8vec2>::value);
static_assert(std::is_trivially_copy_assignable<glm::i16vec2>::value);
static_assert(std::is_trivially_copy_assignable<glm::i32vec2>::value);
static_assert(std::is_trivially_copy_assignable<glm::i64vec2>::value);
static_assert(std::is_trivially_copyable<glm::i8vec2>::value);
static_assert(std::is_trivially_copyable<glm::i16vec2>::value);
static_assert(std::is_trivially_copyable<glm::i32vec2>::value);
static_assert(std::is_trivially_copyable<glm::i64vec2>::value);
static_assert(std::is_copy_constructible<glm::i8vec2>::value);
static_assert(std::is_copy_constructible<glm::i16vec2>::value);
static_assert(std::is_copy_constructible<glm::i32vec2>::value);
static_assert(std::is_copy_constructible<glm::i64vec2>::value);
static_assert(glm::i8vec2::length() == 2);
static_assert(glm::i16vec2::length() == 2);
static_assert(glm::i32vec2::length() == 2);
static_assert(glm::i64vec2::length() == 2);
