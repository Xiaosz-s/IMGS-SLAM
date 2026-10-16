#pragma once

#include "../mat3x3.hpp"
#include "../ext/scalar_uint_sized.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_uint3x3_sized extension included")
#endif

namespace glm
{

	typedef mat<3, 3, uint8, defaultp>				u8mat3x3;

	typedef mat<3, 3, uint16, defaultp>				u16mat3x3;

	typedef mat<3, 3, uint32, defaultp>				u32mat3x3;

	typedef mat<3, 3, uint64, defaultp>				u64mat3x3;

	typedef mat<3, 3, uint8, defaultp>				u8mat3;

	typedef mat<3, 3, uint16, defaultp>				u16mat3;

	typedef mat<3, 3, uint32, defaultp>				u32mat3;

	typedef mat<3, 3, uint64, defaultp>				u64mat3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::u8mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::u16mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::u32mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::u64mat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::u8mat3>::value);
static_assert(std::is_trivially_default_constructible<glm::u16mat3>::value);
static_assert(std::is_trivially_default_constructible<glm::u32mat3>::value);
static_assert(std::is_trivially_default_constructible<glm::u64mat3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::u8mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::u16mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::u32mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::u64mat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::u8mat3>::value);
static_assert(std::is_trivially_copy_assignable<glm::u16mat3>::value);
static_assert(std::is_trivially_copy_assignable<glm::u32mat3>::value);
static_assert(std::is_trivially_copy_assignable<glm::u64mat3>::value);
static_assert(std::is_trivially_copyable<glm::u8mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::u16mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::u32mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::u64mat3x3>::value);
static_assert(std::is_trivially_copyable<glm::u8mat3>::value);
static_assert(std::is_trivially_copyable<glm::u16mat3>::value);
static_assert(std::is_trivially_copyable<glm::u32mat3>::value);
static_assert(std::is_trivially_copyable<glm::u64mat3>::value);
static_assert(std::is_copy_constructible<glm::u8mat3x3>::value);
static_assert(std::is_copy_constructible<glm::u16mat3x3>::value);
static_assert(std::is_copy_constructible<glm::u32mat3x3>::value);
static_assert(std::is_copy_constructible<glm::u64mat3x3>::value);
static_assert(std::is_copy_constructible<glm::u8mat3>::value);
static_assert(std::is_copy_constructible<glm::u16mat3>::value);
static_assert(std::is_copy_constructible<glm::u32mat3>::value);
static_assert(std::is_copy_constructible<glm::u64mat3>::value);
static_assert(glm::u8mat3x3::length() == 3);
static_assert(glm::u16mat3x3::length() == 3);
static_assert(glm::u32mat3x3::length() == 3);
static_assert(glm::u64mat3x3::length() == 3);
static_assert(glm::u8mat3::length() == 3);
static_assert(glm::u16mat3::length() == 3);
static_assert(glm::u32mat3::length() == 3);
static_assert(glm::u64mat3::length() == 3);
