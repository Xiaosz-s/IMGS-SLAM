#pragma once

#include "../mat3x3.hpp"
#include "../ext/scalar_int_sized.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_int3x3_sized extension included")
#endif

namespace glm
{

	typedef mat<3, 3, int8, defaultp>				i8mat3x3;

	typedef mat<3, 3, int16, defaultp>				i16mat3x3;

	typedef mat<3, 3, int32, defaultp>				i32mat3x3;

	typedef mat<3, 3, int64, defaultp>				i64mat3x3;

	typedef mat<3, 3, int8, defaultp>				i8mat3;

	typedef mat<3, 3, int16, defaultp>				i16mat3;

	typedef mat<3, 3, int32, defaultp>				i32mat3;

	typedef mat<3, 3, int64, defaultp>				i64mat3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::i8mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::i16mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::i32mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::i64mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::i8mat3>::value);
static_assert(std::is_trivially_default_constructible<glm::i16mat3>::value);
static_assert(std::is_trivially_default_constructible<glm::i32mat3>::value);
static_assert(std::is_trivially_default_constructible<glm::i64mat3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::i8mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i16mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i32mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i64mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i8mat3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i16mat3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i32mat3>::value);
static_assert(std::is_trivially_copy_assignable<glm::i64mat3>::value);
static_assert(std::is_trivially_copyable<glm::i8mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::i16mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::i32mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::i64mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::i8mat3>::value);
static_assert(std::is_trivially_copyable<glm::i16mat3>::value);
static_assert(std::is_trivially_copyable<glm::i32mat3>::value);
static_assert(std::is_trivially_copyable<glm::i64mat3>::value);
static_assert(std::is_copy_constructible<glm::i8mat3x3>::value);
static_assert(std::is_copy_constructible<glm::i16mat3x3>::value);
static_assert(std::is_copy_constructible<glm::i32mat3x3>::value);
static_assert(std::is_copy_constructible<glm::i64mat3x3>::value);
static_assert(std::is_copy_constructible<glm::i8mat3>::value);
static_assert(std::is_copy_constructible<glm::i16mat3>::value);
static_assert(std::is_copy_constructible<glm::i32mat3>::value);
static_assert(std::is_copy_constructible<glm::i64mat3>::value);
static_assert(glm::i8mat3x3::length() == 3);
static_assert(glm::i16mat3x3::length() == 3);
static_assert(glm::i32mat3x3::length() == 3);
static_assert(glm::i64mat3x3::length() == 3);
static_assert(glm::i8mat3::length() == 3);
static_assert(glm::i16mat3::length() == 3);
static_assert(glm::i32mat3::length() == 3);
static_assert(glm::i64mat3::length() == 3);
