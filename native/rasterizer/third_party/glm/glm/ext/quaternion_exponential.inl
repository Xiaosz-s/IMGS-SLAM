#include "scalar_constants.hpp"

namespace glm
{
	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER qua<T, Q> exp(qua<T, Q> const& q)
	{
		vec<3, T, Q> u(q.x, q.y, q.z);
		T const Angle = glm::length(u);
		if (Angle < epsilon<T>())
			return qua<T, Q>();

		vec<3, T, Q> const v(u / Angle);
		return qua<T, Q>(cos(Angle), sin(Angle) * v);
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER qua<T, Q> log(qua<T, Q> const& q)
	{
		vec<3, T, Q> u(q.x, q.y, q.z);
		T Vec3Len = length(u);

		if (Vec3Len < epsilon<T>())
		{
			if(q.w > static_cast<T>(0))
				return qua<T, Q>::wxyz(log(q.w), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
			else if(q.w < static_cast<T>(0))
				return qua<T, Q>::wxyz(log(-q.w), pi<T>(), static_cast<T>(0), static_cast<T>(0));
			else
				return qua<T, Q>::wxyz(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity());
		}
		else
		{
			T t = atan(Vec3Len, T(q.w)) / Vec3Len;
			T QuatLen2 = Vec3Len * Vec3Len + q.w * q.w;
			return qua<T, Q>::wxyz(static_cast<T>(0.5) * log(QuatLen2), t * q.x, t * q.y, t * q.z);
		}
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER qua<T, Q> pow(qua<T, Q> const& x, T y)
	{

		if(y > -epsilon<T>() && y < epsilon<T>())
			return qua<T, Q>::wxyz(1,0,0,0);

		T magnitude = sqrt(x.x * x.x + x.y * x.y + x.z * x.z + x.w *x.w);

		T Angle;
		if(abs(x.w / magnitude) > cos_one_over_two<T>())
		{

			T VectorMagnitude = x.x * x.x + x.y * x.y + x.z * x.z;

			if (VectorMagnitude < (std::numeric_limits<T>::min)()) {
				
				return qua<T, Q>::wxyz(pow(x.w, y), 0, 0, 0);
			}

			Angle = asin(sqrt(VectorMagnitude) / magnitude);
		}
		else
		{
			
			Angle = acos(x.w / magnitude);
		}

		T NewAngle = Angle * y;
		T Div = sin(NewAngle) / sin(Angle);
		T Mag = pow(magnitude, y - static_cast<T>(1));
		return qua<T, Q>::wxyz(cos(NewAngle) * magnitude * Mag, x.x * Div * Mag, x.y * Div * Mag, x.z * Div * Mag);
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER qua<T, Q> sqrt(qua<T, Q> const& x)
	{
		return pow(x, static_cast<T>(0.5));
	}
}

