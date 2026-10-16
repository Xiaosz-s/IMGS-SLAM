#pragma once

#include "../mat4x3.hpp"
#include "../ext/scalar_int_sized.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_int4x3_sized extension included")
#endif

namespace glm
{

	typedef mat<4, 3, int8, defaultp>				i8mat4x3;

	typedef mat<4, 3, int16, defaultp>				i16mat4x3;

	typedef mat<4, 3, int32, defaultp>				i32mat4x3;

	typedef mat<4, 3, int64, defaultp>				i64mat4x3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::i8mat4x3>::value);
static_assert(std::is_trivially_default_constructible<glm::i16mat4x3>::value);
static_assert(std::is_trivially_default_constructible<glm::i32mat4x3>::value);
static_assert(std::is_trivially_default_constructible<glm::i64mat4x3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::i8mat4x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i16mat4x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i32mat4x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i64mat4x3>::value);
static_assert(std::is_trivially_copyable<glm::i8mat4x3>::value);
static_assert(std::is_trivially_copyable<glm::i16mat4x3>::value);
static_assert(std::is_trivially_copyable<glm::i32mat4x3>::value);
static_assert(std::is_trivially_copyable<glm::i64mat4x3>::value);
static_assert(std::is_copy_constructible<glm::i8mat4x3>::value);
static_assert(std::is_copy_constructible<glm::i16mat4x3>::value);
static_assert(std::is_copy_constructible<glm::i32mat4x3>::value);
static_assert(std::is_copy_constructible<glm::i64mat4x3>::value);
static_assert(glm::i8mat4x3::length() == 4);
static_assert(glm::i16mat4x3::length() == 4);
static_assert(glm::i32mat4x3::length() == 4);
static_assert(glm::i64mat4x3::length() == 4);

