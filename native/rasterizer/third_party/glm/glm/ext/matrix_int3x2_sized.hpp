#pragma once

#include "../mat3x2.hpp"
#include "../ext/scalar_int_sized.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_int3x2_sized extension included")
#endif

namespace glm
{

	typedef mat<3, 2, int8, defaultp>				i8mat3x2;

	typedef mat<3, 2, int16, defaultp>				i16mat3x2;

	typedef mat<3, 2, int32, defaultp>				i32mat3x2;

	typedef mat<3, 2, int64, defaultp>				i64mat3x2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::i8mat3x2>::value);
static_assert(std::is_trivially_default_constructible<glm::i16mat3x2>::value);
static_assert(std::is_trivially_default_constructible<glm::i32mat3x2>::value);
static_assert(std::is_trivially_default_constructible<glm::i64mat3x2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::i8mat3x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::i16mat3x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::i32mat3x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::i64mat3x2>::value);
static_assert(std::is_trivially_copyable<glm::i8mat3x2>::value);
static_assert(std::is_trivially_copyable<glm::i16mat3x2>::value);
static_assert(std::is_trivially_copyable<glm::i32mat3x2>::value);
static_assert(std::is_trivially_copyable<glm::i64mat3x2>::value);
static_assert(std::is_copy_constructible<glm::i8mat3x2>::value);
static_assert(std::is_copy_constructible<glm::i16mat3x2>::value);
static_assert(std::is_copy_constructible<glm::i32mat3x2>::value);
static_assert(std::is_copy_constructible<glm::i64mat3x2>::value);
static_assert(glm::i8mat3x2::length() == 3);
static_assert(glm::i16mat3x2::length() == 3);
static_assert(glm::i32mat3x2::length() == 3);
static_assert(glm::i64mat3x2::length() == 3);
