#ifndef EIGEN_IO_H
#define EIGEN_IO_H

namespace Eigen { 

enum { DontAlignCols = 1 };
enum { StreamPrecision = -1,
       FullPrecision = -2 };

namespace internal {
template<typename Derived>
std::ostream & print_matrix(std::ostream & s, const Derived& _m, const IOFormat& fmt);
}

struct IOFormat
{
  
  IOFormat(int _precision = StreamPrecision, int _flags = 0,
    const std::string& _coeffSeparator = " ",
    const std::string& _rowSeparator = "\n", const std::string& _rowPrefix="", const std::string& _rowSuffix="",
    const std::string& _matPrefix="", const std::string& _matSuffix="", const char _fill=' ')
  : matPrefix(_matPrefix), matSuffix(_matSuffix), rowPrefix(_rowPrefix), rowSuffix(_rowSuffix), rowSeparator(_rowSeparator),
    rowSpacer(""), coeffSeparator(_coeffSeparator), fill(_fill), precision(_precision), flags(_flags)
  {

    if((flags & DontAlignCols))
      return;
    int i = int(matSuffix.length())-1;
    while (i>=0 && matSuffix[i]!='\n')
    {
      rowSpacer += ' ';
      i--;
    }
  }
  std::string matPrefix, matSuffix;
  std::string rowPrefix, rowSuffix, rowSeparator, rowSpacer;
  std::string coeffSeparator;
  char fill;
  int precision;
  int flags;
};

template<typename ExpressionType>
class WithFormat
{
  public:

    WithFormat(const ExpressionType& matrix, const IOFormat& format)
      : m_matrix(matrix), m_format(format)
    {}

    friend std::ostream & operator << (std::ostream & s, const WithFormat& wf)
    {
      return internal::print_matrix(s, wf.m_matrix.eval(), wf.m_format);
    }

  protected:
    typename ExpressionType::Nested m_matrix;
    IOFormat m_format;
};

namespace internal {

template<typename Scalar>
struct significant_decimals_impl
{
  static inline int run()
  {
    return NumTraits<Scalar>::digits10();
  }
};

template<typename Derived>
std::ostream & print_matrix(std::ostream & s, const Derived& _m, const IOFormat& fmt)
{
  using internal::is_same;
  using internal::conditional;

  if(_m.size() == 0)
  {
    s << fmt.matPrefix << fmt.matSuffix;
    return s;
  }
  
  typename Derived::Nested m = _m;
  typedef typename Derived::Scalar Scalar;
  typedef typename
      conditional<
          is_same<Scalar, char>::value ||
            is_same<Scalar, unsigned char>::value ||
            is_same<Scalar, numext::int8_t>::value ||
            is_same<Scalar, numext::uint8_t>::value,
          int,
          typename conditional<
              is_same<Scalar, std::complex<char> >::value ||
                is_same<Scalar, std::complex<unsigned char> >::value ||
                is_same<Scalar, std::complex<numext::int8_t> >::value ||
                is_same<Scalar, std::complex<numext::uint8_t> >::value,
              std::complex<int>,
              const Scalar&
            >::type
        >::type PrintType;

  Index width = 0;

  std::streamsize explicit_precision;
  if(fmt.precision == StreamPrecision)
  {
    explicit_precision = 0;
  }
  else if(fmt.precision == FullPrecision)
  {
    if (NumTraits<Scalar>::IsInteger)
    {
      explicit_precision = 0;
    }
    else
    {
      explicit_precision = significant_decimals_impl<Scalar>::run();
    }
  }
  else
  {
    explicit_precision = fmt.precision;
  }

  std::streamsize old_precision = 0;
  if(explicit_precision) old_precision = s.precision(explicit_precision);

  bool align_cols = !(fmt.flags & DontAlignCols);
  if(align_cols)
  {
    
    for(Index j = 0; j < m.cols(); ++j)
      for(Index i = 0; i < m.rows(); ++i)
      {
        std::stringstream sstr;
        sstr.copyfmt(s);
        sstr << static_cast<PrintType>(m.coeff(i,j));
        width = std::max<Index>(width, Index(sstr.str().length()));
      }
  }
  std::streamsize old_width = s.width();
  char old_fill_character = s.fill();
  s << fmt.matPrefix;
  for(Index i = 0; i < m.rows(); ++i)
  {
    if (i)
      s << fmt.rowSpacer;
    s << fmt.rowPrefix;
    if(width) {
      s.fill(fmt.fill);
      s.width(width);
    }
    s << static_cast<PrintType>(m.coeff(i, 0));
    for(Index j = 1; j < m.cols(); ++j)
    {
      s << fmt.coeffSeparator;
      if(width) {
        s.fill(fmt.fill);
        s.width(width);
      }
      s << static_cast<PrintType>(m.coeff(i, j));
    }
    s << fmt.rowSuffix;
    if( i < m.rows() - 1)
      s << fmt.rowSeparator;
  }
  s << fmt.matSuffix;
  if(explicit_precision) s.precision(old_precision);
  if(width) {
    s.fill(old_fill_character);
    s.width(old_width);
  }
  return s;
}

} 

template<typename Derived>
std::ostream & operator <<
(std::ostream & s,
 const DenseBase<Derived> & m)
{
  return internal::print_matrix(s, m.eval(), EIGEN_DEFAULT_IO_FORMAT);
}

} 

#endif 
