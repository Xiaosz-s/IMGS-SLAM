#pragma once

#include "../glm.hpp"
#include "../ext/scalar_relational.hpp"

#ifndef GLM_ENABLE_EXPERIMENTAL
#	error "GLM: GLM_GTX_pca is an experimental extension and may change in the future. Use #define GLM_ENABLE_EXPERIMENTAL before including it, if you really want to use it."
#elif GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_GTX_pca extension included")
#endif

namespace glm {

	template<length_t D, typename T, qualifier Q>
	GLM_INLINE mat<D, D, T, Q> computeCovarianceMatrix(vec<D, T, Q> const* v, size_t n);

	template<length_t D, typename T, qualifier Q>
	GLM_INLINE mat<D, D, T, Q> computeCovarianceMatrix(vec<D, T, Q> const* v, size_t n, vec<D, T, Q> const& c);

	template<length_t D, typename T, qualifier Q, typename I>
	GLM_FUNC_DECL mat<D, D, T, Q> computeCovarianceMatrix(I const& b, I const& e);

	template<length_t D, typename T, qualifier Q, typename I>
	GLM_FUNC_DECL mat<D, D, T, Q> computeCovarianceMatrix(I const& b, I const& e, vec<D, T, Q> const& c);

	template<length_t D, typename T, qualifier Q>
	GLM_FUNC_DECL unsigned int findEigenvaluesSymReal
	(
		mat<D, D, T, Q> const& covarMat,
		vec<D, T, Q>& outEigenvalues,
		mat<D, D, T, Q>& outEigenvectors
	);

	template<typename T, qualifier Q>
	GLM_FUNC_DISCARD_DECL void sortEigenvalues(vec<2, T, Q>& eigenvalues, mat<2, 2, T, Q>& eigenvectors);

	template<typename T, qualifier Q>
	GLM_FUNC_DISCARD_DECL void sortEigenvalues(vec<3, T, Q>& eigenvalues, mat<3, 3, T, Q>& eigenvectors);

	template<typename T, qualifier Q>
	GLM_FUNC_DISCARD_DECL void sortEigenvalues(vec<4, T, Q>& eigenvalues, mat<4, 4, T, Q>& eigenvectors);

}

#include "pca.inl"
