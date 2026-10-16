namespace glm
{
	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<3, T, Q> slerp
	(
		vec<3, T, Q> const& x,
		vec<3, T, Q> const& y,
		T const& a
	)
	{
		
		T CosAlpha = dot(x, y);
		
		T Alpha = acos(CosAlpha);
		
		T SinAlpha = sin(Alpha);
		
		T t1 = sin((static_cast<T>(1) - a) * Alpha) / SinAlpha;
		T t2 = sin(a * Alpha) / SinAlpha;

		return x * t1 + y * t2;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<2, T, Q> rotate
	(
		vec<2, T, Q> const& v,
		T const& angle
	)
	{
		vec<2, T, Q> Result;
		T const Cos(cos(angle));
		T const Sin(sin(angle));

		Result.x = v.x * Cos - v.y * Sin;
		Result.y = v.x * Sin + v.y * Cos;
		return Result;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<3, T, Q> rotate
	(
		vec<3, T, Q> const& v,
		T const& angle,
		vec<3, T, Q> const& normal
	)
	{
		return mat<3, 3, T, Q>(glm::rotate(angle, normal)) * v;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<4, T, Q> rotate
	(
		vec<4, T, Q> const& v,
		T const& angle,
		vec<3, T, Q> const& normal
	)
	{
		return rotate(angle, normal) * v;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<3, T, Q> rotateX
	(
		vec<3, T, Q> const& v,
		T const& angle
	)
	{
		vec<3, T, Q> Result(v);
		T const Cos(cos(angle));
		T const Sin(sin(angle));

		Result.y = v.y * Cos - v.z * Sin;
		Result.z = v.y * Sin + v.z * Cos;
		return Result;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<3, T, Q> rotateY
	(
		vec<3, T, Q> const& v,
		T const& angle
	)
	{
		vec<3, T, Q> Result = v;
		T const Cos(cos(angle));
		T const Sin(sin(angle));

		Result.x =  v.x * Cos + v.z * Sin;
		Result.z = -v.x * Sin + v.z * Cos;
		return Result;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<3, T, Q> rotateZ
	(
		vec<3, T, Q> const& v,
		T const& angle
	)
	{
		vec<3, T, Q> Result = v;
		T const Cos(cos(angle));
		T const Sin(sin(angle));

		Result.x = v.x * Cos - v.y * Sin;
		Result.y = v.x * Sin + v.y * Cos;
		return Result;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<4, T, Q> rotateX
	(
		vec<4, T, Q> const& v,
		T const& angle
	)
	{
		vec<4, T, Q> Result = v;
		T const Cos(cos(angle));
		T const Sin(sin(angle));

		Result.y = v.y * Cos - v.z * Sin;
		Result.z = v.y * Sin + v.z * Cos;
		return Result;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<4, T, Q> rotateY
	(
		vec<4, T, Q> const& v,
		T const& angle
	)
	{
		vec<4, T, Q> Result = v;
		T const Cos(cos(angle));
		T const Sin(sin(angle));

		Result.x =  v.x * Cos + v.z * Sin;
		Result.z = -v.x * Sin + v.z * Cos;
		return Result;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<4, T, Q> rotateZ
	(
		vec<4, T, Q> const& v,
		T const& angle
	)
	{
		vec<4, T, Q> Result = v;
		T const Cos(cos(angle));
		T const Sin(sin(angle));

		Result.x = v.x * Cos - v.y * Sin;
		Result.y = v.x * Sin + v.y * Cos;
		return Result;
	}

	template<typename T, qualifier Q>
	GLM_FUNC_QUALIFIER mat<4, 4, T, Q> orientation
	(
		vec<3, T, Q> const& Normal,
		vec<3, T, Q> const& Up
	)
	{
		if(all(equal(Normal, Up, epsilon<T>())))
			return mat<4, 4, T, Q>(static_cast<T>(1));

		vec<3, T, Q> RotationAxis = cross(Up, Normal);
		T Angle = acos(dot(Normal, Up));

		return rotate(Angle, RotationAxis);
	}
}
