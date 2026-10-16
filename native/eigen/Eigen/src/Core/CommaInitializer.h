#ifndef EIGEN_COMMAINITIALIZER_H
#define EIGEN_COMMAINITIALIZER_H

namespace Eigen { 

template<typename XprType>
struct CommaInitializer
{
  typedef typename XprType::Scalar Scalar;

  EIGEN_DEVICE_FUNC
  inline CommaInitializer(XprType& xpr, const Scalar& s)
    : m_xpr(xpr), m_row(0), m_col(1), m_currentBlockRows(1)
  {
    eigen_assert(m_xpr.rows() > 0 && m_xpr.cols() > 0
      && "Cannot comma-initialize a 0x0 matrix (operator<<)");
    m_xpr.coeffRef(0,0) = s;
  }

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC
  inline CommaInitializer(XprType& xpr, const DenseBase<OtherDerived>& other)
    : m_xpr(xpr), m_row(0), m_col(other.cols()), m_currentBlockRows(other.rows())
  {
    eigen_assert(m_xpr.rows() >= other.rows() && m_xpr.cols() >= other.cols()
      && "Cannot comma-initialize a 0x0 matrix (operator<<)");
    m_xpr.block(0, 0, other.rows(), other.cols()) = other;
  }

  EIGEN_DEVICE_FUNC
  inline CommaInitializer(const CommaInitializer& o)
  : m_xpr(o.m_xpr), m_row(o.m_row), m_col(o.m_col), m_currentBlockRows(o.m_currentBlockRows) {
    
    const_cast<CommaInitializer&>(o).m_row = m_xpr.rows();
    const_cast<CommaInitializer&>(o).m_col = m_xpr.cols();
    const_cast<CommaInitializer&>(o).m_currentBlockRows = 0;
  }

  EIGEN_DEVICE_FUNC
  CommaInitializer& operator,(const Scalar& s)
  {
    if (m_col==m_xpr.cols())
    {
      m_row+=m_currentBlockRows;
      m_col = 0;
      m_currentBlockRows = 1;
      eigen_assert(m_row<m_xpr.rows()
        && "Too many rows passed to comma initializer (operator<<)");
    }
    eigen_assert(m_col<m_xpr.cols()
      && "Too many coefficients passed to comma initializer (operator<<)");
    eigen_assert(m_currentBlockRows==1);
    m_xpr.coeffRef(m_row, m_col++) = s;
    return *this;
  }

  template<typename OtherDerived>
  EIGEN_DEVICE_FUNC
  CommaInitializer& operator,(const DenseBase<OtherDerived>& other)
  {
    if (m_col==m_xpr.cols() && (other.cols()!=0 || other.rows()!=m_currentBlockRows))
    {
      m_row+=m_currentBlockRows;
      m_col = 0;
      m_currentBlockRows = other.rows();
      eigen_assert(m_row+m_currentBlockRows<=m_xpr.rows()
        && "Too many rows passed to comma initializer (operator<<)");
    }
    eigen_assert((m_col + other.cols() <= m_xpr.cols())
      && "Too many coefficients passed to comma initializer (operator<<)");
    eigen_assert(m_currentBlockRows==other.rows());
    m_xpr.template block<OtherDerived::RowsAtCompileTime, OtherDerived::ColsAtCompileTime>
                    (m_row, m_col, other.rows(), other.cols()) = other;
    m_col += other.cols();
    return *this;
  }

  EIGEN_DEVICE_FUNC
  inline ~CommaInitializer()
#if defined VERIFY_RAISES_ASSERT && (!defined EIGEN_NO_ASSERTION_CHECKING) && defined EIGEN_EXCEPTIONS
  EIGEN_EXCEPTION_SPEC(Eigen::eigen_assert_exception)
#endif
  {
    finished();
  }

  EIGEN_DEVICE_FUNC
  inline XprType& finished() {
      eigen_assert(((m_row+m_currentBlockRows) == m_xpr.rows() || m_xpr.cols() == 0)
           && m_col == m_xpr.cols()
           && "Too few coefficients passed to comma initializer (operator<<)");
      return m_xpr;
  }

  XprType& m_xpr;           
  Index m_row;              
  Index m_col;              
  Index m_currentBlockRows; 
};

template<typename Derived>
EIGEN_DEVICE_FUNC inline CommaInitializer<Derived> DenseBase<Derived>::operator<< (const Scalar& s)
{
  return CommaInitializer<Derived>(*static_cast<Derived*>(this), s);
}

template<typename Derived>
template<typename OtherDerived>
EIGEN_DEVICE_FUNC inline CommaInitializer<Derived>
DenseBase<Derived>::operator<<(const DenseBase<OtherDerived>& other)
{
  return CommaInitializer<Derived>(*static_cast<Derived *>(this), other);
}

} 

#endif 
