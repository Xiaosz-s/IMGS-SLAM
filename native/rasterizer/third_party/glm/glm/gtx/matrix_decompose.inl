#include "../gtc/constants.hpp"
#include "../gtc/epsilon.hpp"
#include "../gtx/transform.hpp"

namespace glm{
namespace detail
{

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<3, T, Q> combine(
		vec<3, T, Q> const& a,
		vec<3, T, Q> const& b,
		T ascl, T bscl)
	{
		return (a * ascl) + (b * bscl);
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<3, T, Q> scale(vec<3, T, Q> const& v, T desiredLength)
	{
		return v * desiredLength / length(v);
	}
}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER bool decompose(mat<4, 4, T, Q> const& ModelMatrix, vec<3, T, Q> & Scale, qua<T, Q> & Orientation, vec<3, T, Q> & Translation, vec<3, T, Q> & Skew, vec<4, T, Q> & Perspective)
	{
		mat<4, 4, T, Q> LocalMatrix(ModelMatrix);

		if(epsilonEqual(LocalMatrix[3][3], static_cast<T>(0), epsilon<T>()))
			return false;

		for(length_t i = 0; i < 4; ++i)
		for(length_t j = 0; j < 4; ++j)
			LocalMatrix[i][j] /= LocalMatrix[3][3];

		mat<4, 4, T, Q> PerspectiveMatrix(LocalMatrix);

		for(length_t i = 0; i < 3; i++)
			PerspectiveMatrix[i][3] = static_cast<T>(0);
		PerspectiveMatrix[3][3] = static_cast<T>(1);

		if(epsilonEqual(determinant(PerspectiveMatrix), static_cast<T>(0), epsilon<T>()))
			return false;

		if(
			epsilonNotEqual(LocalMatrix[0][3], static_cast<T>(0), epsilon<T>()) ||
			epsilonNotEqual(LocalMatrix[1][3], static_cast<T>(0), epsilon<T>()) ||
			epsilonNotEqual(LocalMatrix[2][3], static_cast<T>(0), epsilon<T>()))
		{
			
			vec<4, T, Q> RightHandSide;
			RightHandSide[0] = LocalMatrix[0][3];
			RightHandSide[1] = LocalMatrix[1][3];
			RightHandSide[2] = LocalMatrix[2][3];
			RightHandSide[3] = LocalMatrix[3][3];

			mat<4, 4, T, Q> InversePerspectiveMatrix = glm::inverse(PerspectiveMatrix);
			mat<4, 4, T, Q> TransposedInversePerspectiveMatrix = glm::transpose(InversePerspectiveMatrix);

			Perspective = TransposedInversePerspectiveMatrix * RightHandSide;

			LocalMatrix[0][3] = LocalMatrix[1][3] = LocalMatrix[2][3] = static_cast<T>(0);
			LocalMatrix[3][3] = static_cast<T>(1);
		}
		else
		{
			
			Perspective = vec<4, T, Q>(0, 0, 0, 1);
		}

		Translation = vec<3, T, Q>(LocalMatrix[3]);
		LocalMatrix[3] = vec<4, T, Q>(0, 0, 0, LocalMatrix[3].w);

		vec<3, T, Q> Row[3], Pdum3;

		for(length_t i = 0; i < 3; ++i)
		for(length_t j = 0; j < 3; ++j)
			Row[i][j] = LocalMatrix[i][j];

		Scale.x = length(Row[0]);

		Row[0] = detail::scale(Row[0], static_cast<T>(1));

		Skew.z = dot(Row[0], Row[1]);
		Row[1] = detail::combine(Row[1], Row[0], static_cast<T>(1), -Skew.z);

		Scale.y = length(Row[1]);
		Row[1] = detail::scale(Row[1], static_cast<T>(1));
		Skew.z /= Scale.y;

		Skew.y = glm::dot(Row[0], Row[2]);
		Row[2] = detail::combine(Row[2], Row[0], static_cast<T>(1), -Skew.y);
		Skew.x = glm::dot(Row[1], Row[2]);
		Row[2] = detail::combine(Row[2], Row[1], static_cast<T>(1), -Skew.x);

		Scale.z = length(Row[2]);
		Row[2] = detail::scale(Row[2], static_cast<T>(1));
		Skew.y /= Scale.z;
		Skew.x /= Scale.z;

		Pdum3 = cross(Row[1], Row[2]); 
		if(dot(Row[0], Pdum3) < 0)
		{
			for(length_t i = 0; i < 3; i++)
			{
				Scale[i] *= static_cast<T>(-1);
				Row[i] *= static_cast<T>(-1);
			}
		}

		int i, j, k = 0;
		T root, trace = Row[0].x + Row[1].y + Row[2].z;
		if(trace > static_cast<T>(0))
		{
			root = sqrt(trace + static_cast<T>(1.0));
			Orientation.w = static_cast<T>(0.5) * root;
			root = static_cast<T>(0.5) / root;
			Orientation.x = root * (Row[1].z - Row[2].y);
			Orientation.y = root * (Row[2].x - Row[0].z);
			Orientation.z = root * (Row[0].y - Row[1].x);
		} 
		else
		{
			static int Next[3] = {1, 2, 0};
			i = 0;
			if(Row[1].y > Row[0].x) i = 1;
			if(Row[2].z > Row[i][i]) i = 2;
			j = Next[i];
			k = Next[j];

#           ifdef GLM_FORCE_QUAT_DATA_WXYZ
                int off = 1;
#           else
                int off = 0;
#           endif

			root = sqrt(Row[i][i] - Row[j][j] - Row[k][k] + static_cast<T>(1.0));

			Orientation[i + off] = static_cast<T>(0.5) * root;
			root = static_cast<T>(0.5) / root;
			Orientation[j + off] = root * (Row[i][j] + Row[j][i]);
			Orientation[k + off] = root * (Row[i][k] + Row[k][i]);
			Orientation.w = root * (Row[j][k] - Row[k][j]);
		} 

		return true;
	}

	template <typename T, qualifier Q>
	GLM_FUNC_DECL mat<4, 4, T, Q> recompose(
		vec<3, T, Q> const& scale, qua<T, Q> const& orientation, vec<3, T, Q> const& translation,
		vec<3, T, Q> const& skew, vec<4, T, Q> const& perspective)
	{
		glm::mat4 m = glm::mat4(1.f);

		m[0][3] = perspective.x;
		m[1][3] = perspective.y;
		m[2][3] = perspective.z;
		m[3][3] = perspective.w;

		m *= glm::translate(translation);
		m *= glm::mat4_cast(orientation);

		if (abs(skew.x) > static_cast<T>(0)) {
			glm::mat4 tmp(1.f);
			tmp[2][1] = skew.x;
			m *= tmp;
		}

		if (abs(skew.y) > static_cast<T>(0)) {
			glm::mat4 tmp(1.f);
			tmp[2][0] = skew.y;
			m *= tmp;
		}

		if (abs(skew.z) > static_cast<T>(0)) {
			glm::mat4 tmp(1.f);
			tmp[1][0] = skew.z;
			m *= tmp;
		}

		m *= glm::scale(scale);

		return m;
	}
}
