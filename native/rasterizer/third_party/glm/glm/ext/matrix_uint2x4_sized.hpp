#pragma once

#include "../mat2x4.hpp"
#include "../ext/scalar_uint_sized.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_uint2x4_sized extension included")
#endif

namespace glm
{

	typedef mat<2, 4, uint8, defaultp>				u8mat2x4;

	typedef mat<2, 4, uint16, defaultp>				u16mat2x4;

	typedef mat<2, 4, uint32, defaultp>				u32mat2x4;

	typedef mat<2, 4, uint64, defaultp>				u64mat2x4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::u8mat2x4>::value);
static_assert(std::is_trivially_default_constructible<glm::u16mat2x4>::value);
static_assert(std::is_trivially_default_constructible<glm::u32mat2x4>::value);
static_assert(std::is_trivially_default_constructible<glm::u64mat2x4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::u8mat2x4>::value);
static_assert(std::is_trivially_copy_assignable<glm::u16mat2x4>::value);
static_assert(std::is_trivially_copy_assignable<glm::u32mat2x4>::value);
static_assert(std::is_trivially_copy_assignable<glm::u64mat2x4>::value);
static_assert(std::is_trivially_copyable<glm::u8mat2x4>::value);
static_assert(std::is_trivially_copyable<glm::u16mat2x4>::value);
static_assert(std::is_trivially_copyable<glm::u32mat2x4>::value);
static_assert(std::is_trivially_copyable<glm::u64mat2x4>::value);
static_assert(std::is_copy_constructible<glm::u8mat2x4>::value);
static_assert(std::is_copy_constructible<glm::u16mat2x4>::value);
static_assert(std::is_copy_constructible<glm::u32mat2x4>::value);
static_assert(std::is_copy_constructible<glm::u64mat2x4>::value);
static_assert(glm::u8mat2x4::length() == 2);
static_assert(glm::u16mat2x4::length() == 2);
static_assert(glm::u32mat2x4::length() == 2);
static_assert(glm::u64mat2x4::length() == 2);
