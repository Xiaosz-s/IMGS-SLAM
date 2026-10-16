#ifndef EIGEN_BASIC_PRECONDITIONERS_H
#define EIGEN_BASIC_PRECONDITIONERS_H

namespace Eigen {

template <typename Scalar_>
class DiagonalPreconditioner
{
    typedef Scalar_ Scalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
  public:
    typedef typename Vector::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };

    DiagonalPreconditioner() : m_isInitialized(false) {}

    template<typename MatType>
    explicit DiagonalPreconditioner(const MatType& mat) : m_invdiag(mat.cols())
    {
      compute(mat);
    }

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_invdiag.size(); }
    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_invdiag.size(); }

    template<typename MatType>
    DiagonalPreconditioner& analyzePattern(const MatType& )
    {
      return *this;
    }

    template<typename MatType>
    DiagonalPreconditioner& factorize(const MatType& mat)
    {
      m_invdiag.resize(mat.cols());
      for(int j=0; j<mat.outerSize(); ++j)
      {
        typename MatType::InnerIterator it(mat,j);
        while(it && it.index()!=j) ++it;
        if(it && it.index()==j && it.value()!=Scalar(0))
          m_invdiag(j) = Scalar(1)/it.value();
        else
          m_invdiag(j) = Scalar(1);
      }
      m_isInitialized = true;
      return *this;
    }

    template<typename MatType>
    DiagonalPreconditioner& compute(const MatType& mat)
    {
      return factorize(mat);
    }

    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      x = m_invdiag.array() * b.array() ;
    }

    template<typename Rhs> inline const Solve<DiagonalPreconditioner, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "DiagonalPreconditioner is not initialized.");
      eigen_assert(m_invdiag.size()==b.rows()
                && "DiagonalPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
      return Solve<DiagonalPreconditioner, Rhs>(*this, b.derived());
    }

    ComputationInfo info() { return Success; }

  protected:
    Vector m_invdiag;
    bool m_isInitialized;
};

template <typename Scalar_>
class LeastSquareDiagonalPreconditioner : public DiagonalPreconditioner<Scalar_>
{
    typedef Scalar_ Scalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef DiagonalPreconditioner<Scalar_> Base;
    using Base::m_invdiag;
  public:

    LeastSquareDiagonalPreconditioner() : Base() {}

    template<typename MatType>
    explicit LeastSquareDiagonalPreconditioner(const MatType& mat) : Base()
    {
      compute(mat);
    }

    template<typename MatType>
    LeastSquareDiagonalPreconditioner& analyzePattern(const MatType& )
    {
      return *this;
    }

    template<typename MatType>
    LeastSquareDiagonalPreconditioner& factorize(const MatType& mat)
    {
      
      m_invdiag.resize(mat.cols());
      if(MatType::IsRowMajor)
      {
        m_invdiag.setZero();
        for(Index j=0; j<mat.outerSize(); ++j)
        {
          for(typename MatType::InnerIterator it(mat,j); it; ++it)
            m_invdiag(it.index()) += numext::abs2(it.value());
        }
        for(Index j=0; j<mat.cols(); ++j)
          if(numext::real(m_invdiag(j))>RealScalar(0))
            m_invdiag(j) = RealScalar(1)/numext::real(m_invdiag(j));
      }
      else
      {
        for(Index j=0; j<mat.outerSize(); ++j)
        {
          RealScalar sum = mat.col(j).squaredNorm();
          if(sum>RealScalar(0))
            m_invdiag(j) = RealScalar(1)/sum;
          else
            m_invdiag(j) = RealScalar(1);
        }
      }
      Base::m_isInitialized = true;
      return *this;
    }

    template<typename MatType>
    LeastSquareDiagonalPreconditioner& compute(const MatType& mat)
    {
      return factorize(mat);
    }

    ComputationInfo info() { return Success; }

  protected:
};

class IdentityPreconditioner
{
  public:

    IdentityPreconditioner() {}

    template<typename MatrixType>
    explicit IdentityPreconditioner(const MatrixType& ) {}

    template<typename MatrixType>
    IdentityPreconditioner& analyzePattern(const MatrixType& ) { return *this; }

    template<typename MatrixType>
    IdentityPreconditioner& factorize(const MatrixType& ) { return *this; }

    template<typename MatrixType>
    IdentityPreconditioner& compute(const MatrixType& ) { return *this; }

    template<typename Rhs>
    inline const Rhs& solve(const Rhs& b) const { return b; }

    ComputationInfo info() { return Success; }
};

} 

#endif 
