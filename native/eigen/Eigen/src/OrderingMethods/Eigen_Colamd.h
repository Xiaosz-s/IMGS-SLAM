#ifndef EIGEN_COLAMD_H
#define EIGEN_COLAMD_H

namespace internal {

namespace Colamd {

#ifndef COLAMD_NDEBUG
#define COLAMD_NDEBUG
#endif 

const int NKnobs = 20;

const int NStats = 20;

enum KnobsStatsIndex {
  
  DenseRow = 0,

  DenseCol = 1,

  DefragCount = 2,

  Status = 3,

  Info1 = 4,
  Info2 = 5,
  Info3 = 6
};

enum Status {
  Ok = 0,
  OkButJumbled = 1,
  ErrorANotPresent = -1,
  ErrorPNotPresent = -2,
  ErrorNrowNegative = -3,
  ErrorNcolNegative = -4,
  ErrorNnzNegative = -5,
  ErrorP0Nonzero = -6,
  ErrorATooSmall = -7,
  ErrorColLengthNegative = -8,
  ErrorRowIndexOutOfBounds = -9,
  ErrorOutOfMemory = -10,
  ErrorInternalError = -999
};

template <typename IndexType>
IndexType ones_complement(const IndexType r) {
  return (-(r)-1);
}

const int Empty = -1;

enum RowColumnStatus {
  Alive = 0,
  Dead = -1
};

enum ColumnStatus {
  DeadPrincipal = -1,
  DeadNonPrincipal = -2
};

template <typename IndexType>
struct ColStructure
{
  IndexType start ;   
  
  IndexType length ;  
  union
  {
    IndexType thickness ; 
    
    IndexType parent ;  
    
  } shared1 ;
  union
  {
    IndexType score ; 
    IndexType order ; 
  } shared2 ;
  union
  {
    IndexType headhash ;  
    
    IndexType hash ;  
    IndexType prev ;  
    
  } shared3 ;
  union
  {
    IndexType degree_next ; 
    IndexType hash_next ;   
  } shared4 ;

  inline bool is_dead() const { return start < Alive; }

  inline bool is_alive() const { return start >= Alive; }

  inline bool is_dead_principal() const { return start == DeadPrincipal; }

  inline void kill_principal() { start = DeadPrincipal; }

  inline void kill_non_principal() { start = DeadNonPrincipal; }

};

template <typename IndexType>
struct RowStructure
{
  IndexType start ;   
  IndexType length ;  
  union
  {
    IndexType degree ;  
    IndexType p ;   
  } shared1 ;
  union
  {
    IndexType mark ;  
    IndexType first_column ;
  } shared2 ;

  inline bool is_dead() const { return shared2.mark < Alive; }

  inline bool is_alive() const { return shared2.mark >= Alive; }

  inline void kill() { shared2.mark = Dead; }

};

template <typename IndexType>
inline IndexType colamd_c(IndexType n_col)
{ return IndexType( ((n_col) + 1) * sizeof (ColStructure<IndexType>) / sizeof (IndexType) ) ; }

template <typename IndexType>
inline IndexType  colamd_r(IndexType n_row)
{ return IndexType(((n_row) + 1) * sizeof (RowStructure<IndexType>) / sizeof (IndexType)); }

template <typename IndexType>
static IndexType init_rows_cols (IndexType n_row, IndexType n_col, RowStructure<IndexType> Row [], ColStructure<IndexType> col [], IndexType A [], IndexType p [], IndexType stats[NStats] );

template <typename IndexType>
static void init_scoring (IndexType n_row, IndexType n_col, RowStructure<IndexType> Row [], ColStructure<IndexType> Col [], IndexType A [], IndexType head [], double knobs[NKnobs], IndexType *p_n_row2, IndexType *p_n_col2, IndexType *p_max_deg);

template <typename IndexType>
static IndexType find_ordering (IndexType n_row, IndexType n_col, IndexType Alen, RowStructure<IndexType> Row [], ColStructure<IndexType> Col [], IndexType A [], IndexType head [], IndexType n_col2, IndexType max_deg, IndexType pfree);

template <typename IndexType>
static void order_children (IndexType n_col, ColStructure<IndexType> Col [], IndexType p []);

template <typename IndexType>
static void detect_super_cols (ColStructure<IndexType> Col [], IndexType A [], IndexType head [], IndexType row_start, IndexType row_length ) ;

template <typename IndexType>
static IndexType garbage_collection (IndexType n_row, IndexType n_col, RowStructure<IndexType> Row [], ColStructure<IndexType> Col [], IndexType A [], IndexType *pfree) ;

template <typename IndexType>
static inline  IndexType clear_mark (IndexType n_row, RowStructure<IndexType> Row [] ) ;

#define COLAMD_DEBUG0(params) ;
#define COLAMD_DEBUG1(params) ;
#define COLAMD_DEBUG2(params) ;
#define COLAMD_DEBUG3(params) ;
#define COLAMD_DEBUG4(params) ;

#define COLAMD_ASSERT(expression) ((void) 0)

template <typename IndexType>
inline IndexType recommended ( IndexType nnz, IndexType n_row, IndexType n_col)
{
  if ((nnz) < 0 || (n_row) < 0 || (n_col) < 0)
    return (-1);
  else
    return (2 * (nnz) + colamd_c (n_col) + colamd_r (n_row) + (n_col) + ((nnz) / 5));
}

static inline void set_defaults(double knobs[NKnobs])
{

  int i ;

  if (!knobs)
  {
    return ;      
  }
  for (i = 0 ; i < NKnobs ; i++)
  {
    knobs [i] = 0 ;
  }
  knobs [Colamd::DenseRow] = 0.5 ;  
  knobs [Colamd::DenseCol] = 0.5 ;  
}

template <typename IndexType>
static bool compute_ordering(IndexType n_row, IndexType n_col, IndexType Alen, IndexType *A, IndexType *p, double knobs[NKnobs], IndexType stats[NStats])
{

  IndexType i ;     
  IndexType nnz ;     
  IndexType Row_size ;    
  IndexType Col_size ;    
  IndexType need ;      
  Colamd::RowStructure<IndexType> *Row ;   
  Colamd::ColStructure<IndexType> *Col ;   
  IndexType n_col2 ;    
  IndexType n_row2 ;    
  IndexType ngarbage ;    
  IndexType max_deg ;   
  double default_knobs [NKnobs] ; 

  if (!stats)
  {
    COLAMD_DEBUG0 (("colamd: stats not present\n")) ;
    return (false) ;
  }
  for (i = 0 ; i < NStats ; i++)
  {
    stats [i] = 0 ;
  }
  stats [Colamd::Status] = Colamd::Ok ;
  stats [Colamd::Info1] = -1 ;
  stats [Colamd::Info2] = -1 ;

  if (!A)   
  {
    stats [Colamd::Status] = Colamd::ErrorANotPresent ;
    COLAMD_DEBUG0 (("colamd: A not present\n")) ;
    return (false) ;
  }

  if (!p)   
  {
    stats [Colamd::Status] = Colamd::ErrorPNotPresent ;
    COLAMD_DEBUG0 (("colamd: p not present\n")) ;
    return (false) ;
  }

  if (n_row < 0)  
  {
    stats [Colamd::Status] = Colamd::ErrorNrowNegative ;
    stats [Colamd::Info1] = n_row ;
    COLAMD_DEBUG0 (("colamd: nrow negative %d\n", n_row)) ;
    return (false) ;
  }

  if (n_col < 0)  
  {
    stats [Colamd::Status] = Colamd::ErrorNcolNegative ;
    stats [Colamd::Info1] = n_col ;
    COLAMD_DEBUG0 (("colamd: ncol negative %d\n", n_col)) ;
    return (false) ;
  }

  nnz = p [n_col] ;
  if (nnz < 0)  
  {
    stats [Colamd::Status] = Colamd::ErrorNnzNegative ;
    stats [Colamd::Info1] = nnz ;
    COLAMD_DEBUG0 (("colamd: number of entries negative %d\n", nnz)) ;
    return (false) ;
  }

  if (p [0] != 0)
  {
    stats [Colamd::Status] = Colamd::ErrorP0Nonzero ;
    stats [Colamd::Info1] = p [0] ;
    COLAMD_DEBUG0 (("colamd: p[0] not zero %d\n", p [0])) ;
    return (false) ;
  }

  if (!knobs)
  {
    set_defaults (default_knobs) ;
    knobs = default_knobs ;
  }

  Col_size = colamd_c (n_col) ;
  Row_size = colamd_r (n_row) ;
  need = 2*nnz + n_col + Col_size + Row_size ;

  if (need > Alen)
  {
    
    stats [Colamd::Status] = Colamd::ErrorATooSmall ;
    stats [Colamd::Info1] = need ;
    stats [Colamd::Info2] = Alen ;
    COLAMD_DEBUG0 (("colamd: Need Alen >= %d, given only Alen = %d\n", need,Alen));
    return (false) ;
  }

  Alen -= Col_size + Row_size ;
  Col = (ColStructure<IndexType> *) &A [Alen] ;
  Row = (RowStructure<IndexType> *) &A [Alen + Col_size] ;

  if (!Colamd::init_rows_cols (n_row, n_col, Row, Col, A, p, stats))
  {
    
    COLAMD_DEBUG0 (("colamd: Matrix invalid\n")) ;
    return (false) ;
  }

  Colamd::init_scoring (n_row, n_col, Row, Col, A, p, knobs,
		&n_row2, &n_col2, &max_deg) ;

  ngarbage = Colamd::find_ordering (n_row, n_col, Alen, Row, Col, A, p,
			    n_col2, max_deg, 2*nnz) ;

  Colamd::order_children (n_col, Col, p) ;

  stats [Colamd::DenseRow] = n_row - n_row2 ;
  stats [Colamd::DenseCol] = n_col - n_col2 ;
  stats [Colamd::DefragCount] = ngarbage ;
  COLAMD_DEBUG0 (("colamd: done.\n")) ;
  return (true) ;
}

template <typename IndexType>
static IndexType init_rows_cols  
  (

    IndexType n_row,      
    IndexType n_col,      
    RowStructure<IndexType> Row [],    
    ColStructure<IndexType> Col [],    
    IndexType A [],     
    IndexType p [],     
    IndexType stats [NStats]  
    )
{

  IndexType col ;     
  IndexType row ;     
  IndexType *cp ;     
  IndexType *cp_end ;   
  IndexType *rp ;     
  IndexType *rp_end ;   
  IndexType last_row ;    

  for (col = 0 ; col < n_col ; col++)
  {
    Col [col].start = p [col] ;
    Col [col].length = p [col+1] - p [col] ;

    if ((Col [col].length) < 0) 
    {
      
      stats [Colamd::Status] = Colamd::ErrorColLengthNegative ;
      stats [Colamd::Info1] = col ;
      stats [Colamd::Info2] = Col [col].length ;
      COLAMD_DEBUG0 (("colamd: col %d length %d < 0\n", col, Col [col].length)) ;
      return (false) ;
    }

    Col [col].shared1.thickness = 1 ;
    Col [col].shared2.score = 0 ;
    Col [col].shared3.prev = Empty ;
    Col [col].shared4.degree_next = Empty ;
  }

  stats [Info3] = 0 ;  

  for (row = 0 ; row < n_row ; row++)
  {
    Row [row].length = 0 ;
    Row [row].shared2.mark = -1 ;
  }

  for (col = 0 ; col < n_col ; col++)
  {
    last_row = -1 ;

    cp = &A [p [col]] ;
    cp_end = &A [p [col+1]] ;

    while (cp < cp_end)
    {
      row = *cp++ ;

      if (row < 0 || row >= n_row)
      {
	stats [Colamd::Status] = Colamd::ErrorRowIndexOutOfBounds ;
	stats [Colamd::Info1] = col ;
	stats [Colamd::Info2] = row ;
	stats [Colamd::Info3] = n_row ;
	COLAMD_DEBUG0 (("colamd: row %d col %d out of bounds\n", row, col)) ;
	return (false) ;
      }

      if (row <= last_row || Row [row].shared2.mark == col)
      {

	stats [Colamd::Status] = Colamd::OkButJumbled ;
	stats [Colamd::Info1] = col ;
	stats [Colamd::Info2] = row ;
	(stats [Colamd::Info3]) ++ ;
	COLAMD_DEBUG1 (("colamd: row %d col %d unsorted/duplicate\n",row,col));
      }

      if (Row [row].shared2.mark != col)
      {
	Row [row].length++ ;
      }
      else
      {

	Col [col].length-- ;
      }

      Row [row].shared2.mark = col ;

      last_row = row ;
    }
  }

  Row [0].start = p [n_col] ;
  Row [0].shared1.p = Row [0].start ;
  Row [0].shared2.mark = -1 ;
  for (row = 1 ; row < n_row ; row++)
  {
    Row [row].start = Row [row-1].start + Row [row-1].length ;
    Row [row].shared1.p = Row [row].start ;
    Row [row].shared2.mark = -1 ;
  }

  if (stats [Status] == OkButJumbled)
  {
    
    for (col = 0 ; col < n_col ; col++)
    {
      cp = &A [p [col]] ;
      cp_end = &A [p [col+1]] ;
      while (cp < cp_end)
      {
	row = *cp++ ;
	if (Row [row].shared2.mark != col)
	{
	  A [(Row [row].shared1.p)++] = col ;
	  Row [row].shared2.mark = col ;
	}
      }
    }
  }
  else
  {
    
    for (col = 0 ; col < n_col ; col++)
    {
      cp = &A [p [col]] ;
      cp_end = &A [p [col+1]] ;
      while (cp < cp_end)
      {
	A [(Row [*cp++].shared1.p)++] = col ;
      }
    }
  }

  for (row = 0 ; row < n_row ; row++)
  {
    Row [row].shared2.mark = 0 ;
    Row [row].shared1.degree = Row [row].length ;
  }

  if (stats [Status] == OkButJumbled)
  {
    COLAMD_DEBUG0 (("colamd: reconstructing column form, matrix jumbled\n")) ;

    Col [0].start = 0 ;
    p [0] = Col [0].start ;
    for (col = 1 ; col < n_col ; col++)
    {

      Col [col].start = Col [col-1].start + Col [col-1].length ;
      p [col] = Col [col].start ;
    }

    for (row = 0 ; row < n_row ; row++)
    {
      rp = &A [Row [row].start] ;
      rp_end = rp + Row [row].length ;
      while (rp < rp_end)
      {
	A [(p [*rp++])++] = row ;
      }
    }
  }

  return (true) ;
}

template <typename IndexType>
static void init_scoring
  (

    IndexType n_row,      
    IndexType n_col,      
    RowStructure<IndexType> Row [],    
    ColStructure<IndexType> Col [],    
    IndexType A [],     
    IndexType head [],    
    double knobs [NKnobs],
    IndexType *p_n_row2,    
    IndexType *p_n_col2,    
    IndexType *p_max_deg    
    )
{

  IndexType c ;     
  IndexType r, row ;    
  IndexType *cp ;     
  IndexType deg ;     
  IndexType *cp_end ;   
  IndexType *new_cp ;   
  IndexType col_length ;    
  IndexType score ;     
  IndexType n_col2 ;    
  IndexType n_row2 ;    
  IndexType dense_row_count ; 
  IndexType dense_col_count ; 
  IndexType min_score ;   
  IndexType max_deg ;   
  IndexType next_col ;    

  dense_row_count = numext::maxi(IndexType(0), numext::mini(IndexType(knobs [Colamd::DenseRow] * n_col), n_col)) ;
  dense_col_count = numext::maxi(IndexType(0), numext::mini(IndexType(knobs [Colamd::DenseCol] * n_row), n_row)) ;
  COLAMD_DEBUG1 (("colamd: densecount: %d %d\n", dense_row_count, dense_col_count)) ;
  max_deg = 0 ;
  n_col2 = n_col ;
  n_row2 = n_row ;

  for (c = n_col-1 ; c >= 0 ; c--)
  {
    deg = Col [c].length ;
    if (deg == 0)
    {
      
      Col [c].shared2.order = --n_col2 ;
      Col[c].kill_principal() ;
    }
  }
  COLAMD_DEBUG1 (("colamd: null columns killed: %d\n", n_col - n_col2)) ;

  for (c = n_col-1 ; c >= 0 ; c--)
  {
    
    if (Col[c].is_dead())
    {
      continue ;
    }
    deg = Col [c].length ;
    if (deg > dense_col_count)
    {
      
      Col [c].shared2.order = --n_col2 ;
      
      cp = &A [Col [c].start] ;
      cp_end = cp + Col [c].length ;
      while (cp < cp_end)
      {
	Row [*cp++].shared1.degree-- ;
      }
      Col[c].kill_principal() ;
    }
  }
  COLAMD_DEBUG1 (("colamd: Dense and null columns killed: %d\n", n_col - n_col2)) ;

  for (r = 0 ; r < n_row ; r++)
  {
    deg = Row [r].shared1.degree ;
    COLAMD_ASSERT (deg >= 0 && deg <= n_col) ;
    if (deg > dense_row_count || deg == 0)
    {
      
      Row[r].kill() ;
      --n_row2 ;
    }
    else
    {
      
      max_deg = numext::maxi(max_deg, deg) ;
    }
  }
  COLAMD_DEBUG1 (("colamd: Dense and null rows killed: %d\n", n_row - n_row2)) ;

  for (c = n_col-1 ; c >= 0 ; c--)
  {
    
    if (Col[c].is_dead())
    {
      continue ;
    }
    score = 0 ;
    cp = &A [Col [c].start] ;
    new_cp = cp ;
    cp_end = cp + Col [c].length ;
    while (cp < cp_end)
    {
      
      row = *cp++ ;
      
      if (Row[row].is_dead())
      {
	continue ;
      }
      
      *new_cp++ = row ;
      
      score += Row [row].shared1.degree - 1 ;
      
      score = numext::mini(score, n_col) ;
    }
    
    col_length = (IndexType) (new_cp - &A [Col [c].start]) ;
    if (col_length == 0)
    {

      COLAMD_DEBUG2 (("Newly null killed: %d\n", c)) ;
      Col [c].shared2.order = --n_col2 ;
      Col[c].kill_principal() ;
    }
    else
    {
      
      COLAMD_ASSERT (score >= 0) ;
      COLAMD_ASSERT (score <= n_col) ;
      Col [c].length = col_length ;
      Col [c].shared2.score = score ;
    }
  }
  COLAMD_DEBUG1 (("colamd: Dense, null, and newly-null columns killed: %d\n",
		  n_col-n_col2)) ;

  for (c = 0 ; c <= n_col ; c++)
  {
    head [c] = Empty ;
  }
  min_score = n_col ;

  for (c = n_col-1 ; c >= 0 ; c--)
  {
    
    if (Col[c].is_alive())
    {
      COLAMD_DEBUG4 (("place %d score %d minscore %d ncol %d\n",
		      c, Col [c].shared2.score, min_score, n_col)) ;

      score = Col [c].shared2.score ;

      COLAMD_ASSERT (min_score >= 0) ;
      COLAMD_ASSERT (min_score <= n_col) ;
      COLAMD_ASSERT (score >= 0) ;
      COLAMD_ASSERT (score <= n_col) ;
      COLAMD_ASSERT (head [score] >= Empty) ;

      next_col = head [score] ;
      Col [c].shared3.prev = Empty ;
      Col [c].shared4.degree_next = next_col ;

      if (next_col != Empty)
      {
	Col [next_col].shared3.prev = c ;
      }
      head [score] = c ;

      min_score = numext::mini(min_score, score) ;

    }
  }

  *p_n_col2 = n_col2 ;
  *p_n_row2 = n_row2 ;
  *p_max_deg = max_deg ;
}

template <typename IndexType>
static IndexType find_ordering 
  (

    IndexType n_row,      
    IndexType n_col,      
    IndexType Alen,     
    RowStructure<IndexType> Row [],    
    ColStructure<IndexType> Col [],    
    IndexType A [],     
    IndexType head [],    
    IndexType n_col2,     
    IndexType max_deg,    
    IndexType pfree     
    )
{

  IndexType k ;     
  IndexType pivot_col ;   
  IndexType *cp ;     
  IndexType *rp ;     
  IndexType pivot_row ;   
  IndexType *new_cp ;   
  IndexType *new_rp ;   
  IndexType pivot_row_start ; 
  IndexType pivot_row_degree ;  
  IndexType pivot_row_length ;  
  IndexType pivot_col_score ; 
  IndexType needed_memory ;   
  IndexType *cp_end ;   
  IndexType *rp_end ;   
  IndexType row ;     
  IndexType col ;     
  IndexType max_score ;   
  IndexType cur_score ;   
  unsigned int hash ;   
  IndexType head_column ;   
  IndexType first_col ;   
  IndexType tag_mark ;    
  IndexType row_mark ;    
  IndexType set_difference ;  
  IndexType min_score ;   
  IndexType col_thickness ;   
  IndexType max_mark ;    
  IndexType pivot_col_thickness ; 
  IndexType prev_col ;    
  IndexType next_col ;    
  IndexType ngarbage ;    

  max_mark = INT_MAX - n_col ;  
  tag_mark = Colamd::clear_mark (n_row, Row) ;
  min_score = 0 ;
  ngarbage = 0 ;
  COLAMD_DEBUG1 (("colamd: Ordering, n_col2=%d\n", n_col2)) ;

  for (k = 0 ; k < n_col2 ; )
  {

    COLAMD_ASSERT (min_score >= 0) ;
    COLAMD_ASSERT (min_score <= n_col) ;
    COLAMD_ASSERT (head [min_score] >= Empty) ;

    while (min_score < n_col && head [min_score] == Empty)
    {
      min_score++ ;
    }
    pivot_col = head [min_score] ;
    COLAMD_ASSERT (pivot_col >= 0 && pivot_col <= n_col) ;
    next_col = Col [pivot_col].shared4.degree_next ;
    head [min_score] = next_col ;
    if (next_col != Empty)
    {
      Col [next_col].shared3.prev = Empty ;
    }

    COLAMD_ASSERT (Col[pivot_col].is_alive()) ;
    COLAMD_DEBUG3 (("Pivot col: %d\n", pivot_col)) ;

    pivot_col_score = Col [pivot_col].shared2.score ;

    Col [pivot_col].shared2.order = k ;

    pivot_col_thickness = Col [pivot_col].shared1.thickness ;
    k += pivot_col_thickness ;
    COLAMD_ASSERT (pivot_col_thickness > 0) ;

    needed_memory = numext::mini(pivot_col_score, n_col - k) ;
    if (pfree + needed_memory >= Alen)
    {
      pfree = Colamd::garbage_collection (n_row, n_col, Row, Col, A, &A [pfree]) ;
      ngarbage++ ;
      
      COLAMD_ASSERT (pfree + needed_memory < Alen) ;
      
      tag_mark = Colamd::clear_mark (n_row, Row) ;

    }

    pivot_row_start = pfree ;

    pivot_row_degree = 0 ;

    Col [pivot_col].shared1.thickness = -pivot_col_thickness ;

    cp = &A [Col [pivot_col].start] ;
    cp_end = cp + Col [pivot_col].length ;
    while (cp < cp_end)
    {
      
      row = *cp++ ;
      COLAMD_DEBUG4 (("Pivot col pattern %d %d\n", Row[row].is_alive(), row)) ;
      
      if (Row[row].is_dead())
      {
	continue ;
      }
      rp = &A [Row [row].start] ;
      rp_end = rp + Row [row].length ;
      while (rp < rp_end)
      {
	
	col = *rp++ ;
	
	col_thickness = Col [col].shared1.thickness ;
	if (col_thickness > 0 && Col[col].is_alive())
	{
	  
	  Col [col].shared1.thickness = -col_thickness ;
	  COLAMD_ASSERT (pfree < Alen) ;
	  
	  A [pfree++] = col ;
	  pivot_row_degree += col_thickness ;
	}
      }
    }

    Col [pivot_col].shared1.thickness = pivot_col_thickness ;
    max_deg = numext::maxi(max_deg, pivot_row_degree) ;

    cp = &A [Col [pivot_col].start] ;
    cp_end = cp + Col [pivot_col].length ;
    while (cp < cp_end)
    {
      
      row = *cp++ ;
      COLAMD_DEBUG3 (("Kill row in pivot col: %d\n", row)) ;
      Row[row].kill() ;
    }

    pivot_row_length = pfree - pivot_row_start ;
    if (pivot_row_length > 0)
    {
      
      pivot_row = A [Col [pivot_col].start] ;
      COLAMD_DEBUG3 (("Pivotal row is %d\n", pivot_row)) ;
    }
    else
    {
      
      pivot_row = Empty ;
      COLAMD_ASSERT (pivot_row_length == 0) ;
    }
    COLAMD_ASSERT (Col [pivot_col].length > 0 || pivot_row_length == 0) ;

    COLAMD_DEBUG3 (("** Computing set differences phase. **\n")) ;

    COLAMD_DEBUG3 (("Pivot row: ")) ;
    
    rp = &A [pivot_row_start] ;
    rp_end = rp + pivot_row_length ;
    while (rp < rp_end)
    {
      col = *rp++ ;
      COLAMD_ASSERT (Col[col].is_alive() && col != pivot_col) ;
      COLAMD_DEBUG3 (("Col: %d\n", col)) ;

      col_thickness = -Col [col].shared1.thickness ;
      COLAMD_ASSERT (col_thickness > 0) ;
      Col [col].shared1.thickness = col_thickness ;

      cur_score = Col [col].shared2.score ;
      prev_col = Col [col].shared3.prev ;
      next_col = Col [col].shared4.degree_next ;
      COLAMD_ASSERT (cur_score >= 0) ;
      COLAMD_ASSERT (cur_score <= n_col) ;
      COLAMD_ASSERT (cur_score >= Empty) ;
      if (prev_col == Empty)
      {
	head [cur_score] = next_col ;
      }
      else
      {
	Col [prev_col].shared4.degree_next = next_col ;
      }
      if (next_col != Empty)
      {
	Col [next_col].shared3.prev = prev_col ;
      }

      cp = &A [Col [col].start] ;
      cp_end = cp + Col [col].length ;
      while (cp < cp_end)
      {
	
	row = *cp++ ;
	
	if (Row[row].is_dead())
	{
	  continue ;
	}
  row_mark = Row [row].shared2.mark ;
	COLAMD_ASSERT (row != pivot_row) ;
	set_difference = row_mark - tag_mark ;
	
	if (set_difference < 0)
	{
	  COLAMD_ASSERT (Row [row].shared1.degree <= max_deg) ;
	  set_difference = Row [row].shared1.degree ;
	}
	
	set_difference -= col_thickness ;
	COLAMD_ASSERT (set_difference >= 0) ;
	
	if (set_difference == 0)
	{
	  COLAMD_DEBUG3 (("aggressive absorption. Row: %d\n", row)) ;
	  Row[row].kill() ;
	}
	else
	{
	  
	  Row [row].shared2.mark = set_difference + tag_mark ;
	}
      }
    }

    COLAMD_DEBUG3 (("** Adding set differences phase. **\n")) ;

    rp = &A [pivot_row_start] ;
    rp_end = rp + pivot_row_length ;
    while (rp < rp_end)
    {
      
      col = *rp++ ;
      COLAMD_ASSERT (Col[col].is_alive() && col != pivot_col) ;
      hash = 0 ;
      cur_score = 0 ;
      cp = &A [Col [col].start] ;
      
      new_cp = cp ;
      cp_end = cp + Col [col].length ;

      COLAMD_DEBUG4 (("Adding set diffs for Col: %d.\n", col)) ;

      while (cp < cp_end)
      {
	
	row = *cp++ ;
	COLAMD_ASSERT(row >= 0 && row < n_row) ;
	
	if (Row [row].is_dead())
	{
	  continue ;
	}
  row_mark = Row [row].shared2.mark ;
	COLAMD_ASSERT (row_mark > tag_mark) ;
	
	*new_cp++ = row ;
	
	hash += row ;
	
	cur_score += row_mark - tag_mark ;
	
	cur_score = numext::mini(cur_score, n_col) ;
      }

      Col [col].length = (IndexType) (new_cp - &A [Col [col].start]) ;

      if (Col [col].length == 0)
      {
	COLAMD_DEBUG4 (("further mass elimination. Col: %d\n", col)) ;
	
	Col[col].kill_principal() ;
	pivot_row_degree -= Col [col].shared1.thickness ;
	COLAMD_ASSERT (pivot_row_degree >= 0) ;
	
	Col [col].shared2.order = k ;
	
	k += Col [col].shared1.thickness ;
      }
      else
      {

	COLAMD_DEBUG4 (("Preparing supercol detection for Col: %d.\n", col)) ;

	Col [col].shared2.score = cur_score ;

	hash %= n_col + 1 ;

	COLAMD_DEBUG4 ((" Hash = %d, n_col = %d.\n", hash, n_col)) ;
	COLAMD_ASSERT (hash <= n_col) ;

	head_column = head [hash] ;
	if (head_column > Empty)
	{

	  first_col = Col [head_column].shared3.headhash ;
	  Col [head_column].shared3.headhash = col ;
	}
	else
	{
	  
	  first_col = - (head_column + 2) ;
	  head [hash] = - (col + 2) ;
	}
	Col [col].shared4.hash_next = first_col ;

	Col [col].shared3.hash = (IndexType) hash ;
	COLAMD_ASSERT (Col[col].is_alive()) ;
      }
    }

    COLAMD_DEBUG3 (("** Supercolumn detection phase. **\n")) ;

    Colamd::detect_super_cols (Col, A, head, pivot_row_start, pivot_row_length) ;

    Col[pivot_col].kill_principal() ;

    tag_mark += (max_deg + 1) ;
    if (tag_mark >= max_mark)
    {
      COLAMD_DEBUG2 (("clearing tag_mark\n")) ;
      tag_mark = Colamd::clear_mark (n_row, Row) ;
    }

    COLAMD_DEBUG3 (("** Finalize scores phase. **\n")) ;

    rp = &A [pivot_row_start] ;
    
    new_rp = rp ;
    rp_end = rp + pivot_row_length ;
    while (rp < rp_end)
    {
      col = *rp++ ;
      
      if (Col[col].is_dead())
      {
	continue ;
      }
      *new_rp++ = col ;
      
      A [Col [col].start + (Col [col].length++)] = pivot_row ;

      cur_score = Col [col].shared2.score + pivot_row_degree ;

      max_score = n_col - k - Col [col].shared1.thickness ;

      cur_score -= Col [col].shared1.thickness ;

      cur_score = numext::mini(cur_score, max_score) ;
      COLAMD_ASSERT (cur_score >= 0) ;

      Col [col].shared2.score = cur_score ;

      COLAMD_ASSERT (min_score >= 0) ;
      COLAMD_ASSERT (min_score <= n_col) ;
      COLAMD_ASSERT (cur_score >= 0) ;
      COLAMD_ASSERT (cur_score <= n_col) ;
      COLAMD_ASSERT (head [cur_score] >= Empty) ;
      next_col = head [cur_score] ;
      Col [col].shared4.degree_next = next_col ;
      Col [col].shared3.prev = Empty ;
      if (next_col != Empty)
      {
	Col [next_col].shared3.prev = col ;
      }
      head [cur_score] = col ;

      min_score = numext::mini(min_score, cur_score) ;

    }

    if (pivot_row_degree > 0)
    {

      Row [pivot_row].start  = pivot_row_start ;
      Row [pivot_row].length = (IndexType) (new_rp - &A[pivot_row_start]) ;
      Row [pivot_row].shared1.degree = pivot_row_degree ;
      Row [pivot_row].shared2.mark = 0 ;
      
    }
  }

  return (ngarbage) ;
}

template <typename IndexType>
static inline  void order_children
(

  IndexType n_col,      
  ColStructure<IndexType> Col [],    
  IndexType p []      
  )
{

  IndexType i ;     
  IndexType c ;     
  IndexType parent ;    
  IndexType order ;     

  for (i = 0 ; i < n_col ; i++)
  {
    
    COLAMD_ASSERT (col_is_dead(Col, i)) ;
    if (!Col[i].is_dead_principal() && Col [i].shared2.order == Empty)
    {
      parent = i ;
      
      do
      {
	parent = Col [parent].shared1.parent ;
      } while (!Col[parent].is_dead_principal()) ;

      c = i ;
      
      order = Col [parent].shared2.order ;

      do
      {
	COLAMD_ASSERT (Col [c].shared2.order == Empty) ;

	Col [c].shared2.order = order++ ;
	
	Col [c].shared1.parent = parent ;

	c = Col [c].shared1.parent ;

      } while (Col [c].shared2.order == Empty) ;

      Col [parent].shared2.order = order ;
    }
  }

  for (c = 0 ; c < n_col ; c++)
  {
    p [Col [c].shared2.order] = c ;
  }
}

template <typename IndexType>
static void detect_super_cols
(

  ColStructure<IndexType> Col [],    
  IndexType A [],     
  IndexType head [],    
  IndexType row_start,    
  IndexType row_length    
)
{

  IndexType hash ;      
  IndexType *rp ;     
  IndexType c ;     
  IndexType super_c ;   
  IndexType *cp1 ;      
  IndexType *cp2 ;      
  IndexType length ;    
  IndexType prev_c ;    
  IndexType i ;     
  IndexType *rp_end ;   
  IndexType col ;     
  IndexType head_column ;   
  IndexType first_col ;   

  rp = &A [row_start] ;
  rp_end = rp + row_length ;
  while (rp < rp_end)
  {
    col = *rp++ ;
    if (Col[col].is_dead())
    {
      continue ;
    }

    hash = Col [col].shared3.hash ;
    COLAMD_ASSERT (hash <= n_col) ;

    head_column = head [hash] ;
    if (head_column > Empty)
    {
      first_col = Col [head_column].shared3.headhash ;
    }
    else
    {
      first_col = - (head_column + 2) ;
    }

    for (super_c = first_col ; super_c != Empty ;
	 super_c = Col [super_c].shared4.hash_next)
    {
      COLAMD_ASSERT (Col [super_c].is_alive()) ;
      COLAMD_ASSERT (Col [super_c].shared3.hash == hash) ;
      length = Col [super_c].length ;

      prev_c = super_c ;

      for (c = Col [super_c].shared4.hash_next ;
	   c != Empty ; c = Col [c].shared4.hash_next)
      {
	COLAMD_ASSERT (c != super_c) ;
	COLAMD_ASSERT (Col[c].is_alive()) ;
	COLAMD_ASSERT (Col [c].shared3.hash == hash) ;

	if (Col [c].length != length ||
	    Col [c].shared2.score != Col [super_c].shared2.score)
	{
	  prev_c = c ;
	  continue ;
	}

	cp1 = &A [Col [super_c].start] ;
	cp2 = &A [Col [c].start] ;

	for (i = 0 ; i < length ; i++)
	{
	  
	  COLAMD_ASSERT ( cp1->is_alive() );
	  COLAMD_ASSERT ( cp2->is_alive() );

	  if (*cp1++ != *cp2++)
	  {
	    break ;
	  }
	}

	if (i != length)
	{
	  prev_c = c ;
	  continue ;
	}

	COLAMD_ASSERT (Col [c].shared2.score == Col [super_c].shared2.score) ;

	Col [super_c].shared1.thickness += Col [c].shared1.thickness ;
	Col [c].shared1.parent = super_c ;
	Col[c].kill_non_principal() ;
	
	Col [c].shared2.order = Empty ;
	
	Col [prev_c].shared4.hash_next = Col [c].shared4.hash_next ;
      }
    }

    if (head_column > Empty)
    {
      
      Col [head_column].shared3.headhash = Empty ;
    }
    else
    {
      
      head [hash] = Empty ;
    }
  }
}

template <typename IndexType>
static IndexType garbage_collection  
  (

    IndexType n_row,      
    IndexType n_col,      
    RowStructure<IndexType> Row [],    
    ColStructure<IndexType> Col [],    
    IndexType A [],     
    IndexType *pfree      
    )
{

  IndexType *psrc ;     
  IndexType *pdest ;    
  IndexType j ;     
  IndexType r ;     
  IndexType c ;     
  IndexType length ;    

  pdest = &A[0] ;
  for (c = 0 ; c < n_col ; c++)
  {
    if (Col[c].is_alive())
    {
      psrc = &A [Col [c].start] ;

      COLAMD_ASSERT (pdest <= psrc) ;
      Col [c].start = (IndexType) (pdest - &A [0]) ;
      length = Col [c].length ;
      for (j = 0 ; j < length ; j++)
      {
	r = *psrc++ ;
	if (Row[r].is_alive())
	{
	  *pdest++ = r ;
	}
      }
      Col [c].length = (IndexType) (pdest - &A [Col [c].start]) ;
    }
  }

  for (r = 0 ; r < n_row ; r++)
  {
    if (Row[r].is_alive())
    {
      if (Row [r].length == 0)
      {
        
        COLAMD_DEBUG3 (("Defrag row kill\n")) ;
        Row[r].kill() ;
      }
      else
      {
        
        psrc = &A [Row [r].start] ;
        Row [r].shared2.first_column = *psrc ;
        COLAMD_ASSERT (Row[r].is_alive()) ;
        
        *psrc = ones_complement(r) ;

      }
    }
  }

  psrc = pdest ;
  while (psrc < pfree)
  {
    
    if (*psrc++ < 0)
    {
      psrc-- ;
      
      r = ones_complement(*psrc) ;
      COLAMD_ASSERT (r >= 0 && r < n_row) ;
      
      *psrc = Row [r].shared2.first_column ;
      COLAMD_ASSERT (Row[r].is_alive()) ;

      COLAMD_ASSERT (pdest <= psrc) ;
      Row [r].start = (IndexType) (pdest - &A [0]) ;
      length = Row [r].length ;
      for (j = 0 ; j < length ; j++)
      {
	c = *psrc++ ;
	if (Col[c].is_alive())
	{
	  *pdest++ = c ;
	}
      }
      Row [r].length = (IndexType) (pdest - &A [Row [r].start]) ;

    }
  }
  
  COLAMD_ASSERT (debug_rows == 0) ;

  return ((IndexType) (pdest - &A [0])) ;
}

template <typename IndexType>
static inline  IndexType clear_mark  
  (

    IndexType n_row,    
    RowStructure<IndexType> Row [] 
    )
{

  IndexType r ;

  for (r = 0 ; r < n_row ; r++)
  {
    if (Row[r].is_alive())
    {
      Row [r].shared2.mark = 0 ;
    }
  }
  return (1) ;
}

} 

} 
#endif
