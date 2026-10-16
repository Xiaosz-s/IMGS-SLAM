#pragma once

#include "../detail/setup.hpp"
#include "../detail/qualifier.hpp"

#ifndef GLM_ENABLE_EXPERIMENTAL
#	error "GLM: GLM_GTX_exterior_product is an experimental extension and may change in the future. Use #define GLM_ENABLE_EXPERIMENTAL before including it, if you really want to use it."
#elif GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_GTX_exterior_product extension included")
#endif

namespace glm
{

	template<typename T, qualifier Q>
	GLM_FUNC_DECL GLM_CONSTEXPR T cross(vec<2, T, Q> const& v, vec<2, T, Q> const& u);

} 

#include "exterior_product.inl"
