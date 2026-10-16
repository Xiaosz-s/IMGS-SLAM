#if defined(EIGEN_USE_SYCL) && \
    !defined(EIGEN_CXX11_TENSOR_TENSOR_SYCL_STORAGE_MEMORY_H)
#define EIGEN_CXX11_TENSOR_TENSOR_SYCL_STORAGE_MEMORY_H

#include <CL/sycl.hpp>
#ifdef EIGEN_EXCEPTIONS
#include <stdexcept>
#endif
#include <cstddef>
#include <queue>
#include <set>
#include <unordered_map>

namespace Eigen {
namespace TensorSycl {
namespace internal {

using sycl_acc_target = cl::sycl::access::target;
using sycl_acc_mode = cl::sycl::access::mode;

using buffer_data_type_t = uint8_t;
const sycl_acc_target default_acc_target = sycl_acc_target::global_buffer;
const sycl_acc_mode default_acc_mode = sycl_acc_mode::read_write;

class PointerMapper {
 public:
  using base_ptr_t = std::intptr_t;

  struct virtual_pointer_t {

    base_ptr_t m_contents;

    operator void *() const { return reinterpret_cast<void *>(m_contents); }

    operator base_ptr_t() const { return m_contents; }

    virtual_pointer_t operator+(size_t off) { return m_contents + off; }

    bool operator<(virtual_pointer_t rhs) const {
      return (static_cast<base_ptr_t>(m_contents) <
              static_cast<base_ptr_t>(rhs.m_contents));
    }

    bool operator>(virtual_pointer_t rhs) const {
      return (static_cast<base_ptr_t>(m_contents) >
              static_cast<base_ptr_t>(rhs.m_contents));
    }

    bool operator==(virtual_pointer_t rhs) const {
      return (static_cast<base_ptr_t>(m_contents) ==
              static_cast<base_ptr_t>(rhs.m_contents));
    }

    bool operator!=(virtual_pointer_t rhs) const {
      return !(this->operator==(rhs));
    }

    virtual_pointer_t(const void *ptr)
        : m_contents(reinterpret_cast<base_ptr_t>(ptr)){};

    virtual_pointer_t(base_ptr_t u) : m_contents(u){};
  };

  const virtual_pointer_t null_virtual_ptr = nullptr;

  static inline bool is_nullptr(virtual_pointer_t ptr) {
    return (static_cast<void *>(ptr) == nullptr);
  }

  using buffer_t = cl::sycl::buffer_mem;

  struct pMapNode_t {
    buffer_t m_buffer;
    size_t m_size;
    bool m_free;

    pMapNode_t(buffer_t b, size_t size, bool f)
        : m_buffer{b}, m_size{size}, m_free{f} {
      m_buffer.set_final_data(nullptr);
    }

    bool operator<=(const pMapNode_t &rhs) { return (m_size <= rhs.m_size); }
  };

  using pointerMap_t = std::map<virtual_pointer_t, pMapNode_t>;

  typename pointerMap_t::iterator get_insertion_point(size_t requiredSize) {
    typename pointerMap_t::iterator retVal;
    bool reuse = false;
    if (!m_freeList.empty()) {
      
      for (auto freeElem : m_freeList) {
        if (freeElem->second.m_size >= requiredSize) {
          retVal = freeElem;
          reuse = true;
          
          m_freeList.erase(freeElem);
          break;
        }
      }
    }
    if (!reuse) {
      retVal = std::prev(m_pointerMap.end());
    }
    return retVal;
  }

  typename pointerMap_t::iterator get_node(const virtual_pointer_t ptr) {
    if (this->count() == 0) {
      m_pointerMap.clear();
      EIGEN_THROW_X(std::out_of_range("There are no pointers allocated\n"));

    }
    if (is_nullptr(ptr)) {
      m_pointerMap.clear();
      EIGEN_THROW_X(std::out_of_range("Cannot access null pointer\n"));
    }

    auto node = m_pointerMap.lower_bound(ptr);

    if (node == std::end(m_pointerMap)) {
      --node;
    } else if (node->first != ptr) {
      if (node == std::begin(m_pointerMap)) {
        m_pointerMap.clear();
        EIGEN_THROW_X(
            std::out_of_range("The pointer is not registered in the map\n"));

      }
      --node;
    }

    return node;
  }

  template <typename buffer_data_type = buffer_data_type_t>
  cl::sycl::buffer<buffer_data_type, 1> get_buffer(
      const virtual_pointer_t ptr) {
    using sycl_buffer_t = cl::sycl::buffer<buffer_data_type, 1>;

    auto node = get_node(ptr);
    eigen_assert(node->first == ptr || node->first < ptr);
    eigen_assert(ptr < static_cast<virtual_pointer_t>(node->second.m_size +
                                                      node->first));
    return *(static_cast<sycl_buffer_t *>(&node->second.m_buffer));
  }

  template <sycl_acc_mode access_mode = default_acc_mode,
            sycl_acc_target access_target = default_acc_target,
            typename buffer_data_type = buffer_data_type_t>
  cl::sycl::accessor<buffer_data_type, 1, access_mode, access_target>
  get_access(const virtual_pointer_t ptr) {
    auto buf = get_buffer<buffer_data_type>(ptr);
    return buf.template get_access<access_mode, access_target>();
  }

  template <sycl_acc_mode access_mode = default_acc_mode,
            sycl_acc_target access_target = default_acc_target,
            typename buffer_data_type = buffer_data_type_t>
  cl::sycl::accessor<buffer_data_type, 1, access_mode, access_target>
  get_access(const virtual_pointer_t ptr, cl::sycl::handler &cgh) {
    auto buf = get_buffer<buffer_data_type>(ptr);
    return buf.template get_access<access_mode, access_target>(cgh);
  }

  inline std::ptrdiff_t get_offset(const virtual_pointer_t ptr) {

    auto node = get_node(ptr);
    auto start = node->first;
    eigen_assert(start == ptr || start < ptr);
    eigen_assert(ptr < start + node->second.m_size);
    return (ptr - start);
  }

  template <typename buffer_data_type>
  inline size_t get_element_offset(const virtual_pointer_t ptr) {
    return get_offset(ptr) / sizeof(buffer_data_type);
  }

  PointerMapper(base_ptr_t baseAddress = 4096)
      : m_pointerMap{}, m_freeList{}, m_baseAddress{baseAddress} {
    if (m_baseAddress == 0) {
      EIGEN_THROW_X(std::invalid_argument("Base address cannot be zero\n"));
    }
  };

  PointerMapper(const PointerMapper &) = delete;

  inline void clear() {
    m_freeList.clear();
    m_pointerMap.clear();
  }

  inline virtual_pointer_t add_pointer(const buffer_t &b) {
    return add_pointer_impl(b);
  }

  inline virtual_pointer_t add_pointer(buffer_t &&b) {
    return add_pointer_impl(b);
  }

  void fuse_forward(typename pointerMap_t::iterator &node) {
    while (node != std::prev(m_pointerMap.end())) {

      auto fwd_node = std::next(node);
      if (!fwd_node->second.m_free) {
        break;
      }
      auto fwd_size = fwd_node->second.m_size;
      m_freeList.erase(fwd_node);
      m_pointerMap.erase(fwd_node);

      node->second.m_size += fwd_size;
    }
  }

  void fuse_backward(typename pointerMap_t::iterator &node) {
    while (node != m_pointerMap.begin()) {

      auto prev_node = std::prev(node);
      if (!prev_node->second.m_free) {
        break;
      }
      prev_node->second.m_size += node->second.m_size;

      m_freeList.erase(node);
      m_pointerMap.erase(node);

      node = prev_node;
    }
  }

  template <bool ReUse = true>
  void remove_pointer(const virtual_pointer_t ptr) {
    if (is_nullptr(ptr)) {
      return;
    }
    auto node = this->get_node(ptr);

    node->second.m_free = true;
    m_freeList.emplace(node);

    fuse_forward(node);
    fuse_backward(node);

    if (node == std::prev(m_pointerMap.end())) {
      m_freeList.erase(node);
      m_pointerMap.erase(node);
    }
  }

  size_t count() const { return (m_pointerMap.size() - m_freeList.size()); }

 private:

  template <class BufferT>
  virtual_pointer_t add_pointer_impl(BufferT b) {
    virtual_pointer_t retVal = nullptr;
    size_t bufSize = b.get_count();
    pMapNode_t p{b, bufSize, false};
    
    if (m_pointerMap.empty()) {
      virtual_pointer_t initialVal{m_baseAddress};
      m_pointerMap.emplace(initialVal, p);
      return initialVal;
    }

    auto lastElemIter = get_insertion_point(bufSize);
    
    if (lastElemIter->second.m_free) {
      lastElemIter->second.m_buffer = b;
      lastElemIter->second.m_free = false;

      if (lastElemIter->second.m_size > bufSize) {
        
        auto remainingSize = lastElemIter->second.m_size - bufSize;
        pMapNode_t p2{b, remainingSize, true};

        lastElemIter->second.m_size = bufSize;

        auto newFreePtr = lastElemIter->first + bufSize;
        auto freeNode = m_pointerMap.emplace(newFreePtr, p2).first;
        m_freeList.emplace(freeNode);
      }

      retVal = lastElemIter->first;
    } else {
      size_t lastSize = lastElemIter->second.m_size;
      retVal = lastElemIter->first + lastSize;
      m_pointerMap.emplace(retVal, p);
    }
    return retVal;
  }

  struct SortBySize {
    bool operator()(typename pointerMap_t::iterator a,
                    typename pointerMap_t::iterator b) const {
      return ((a->first < b->first) && (a->second <= b->second)) ||
             ((a->first < b->first) && (b->second <= a->second));
    }
  };

  pointerMap_t m_pointerMap;

  std::set<typename pointerMap_t::iterator, SortBySize> m_freeList;

  std::intptr_t m_baseAddress;
};

template <>
inline void PointerMapper::remove_pointer<false>(const virtual_pointer_t ptr) {
  if (is_nullptr(ptr)) {
    return;
  }
  m_pointerMap.erase(this->get_node(ptr));
}

inline void *SYCLmalloc(size_t size, PointerMapper &pMap) {
  if (size == 0) {
    return nullptr;
  }
  
  using buffer_t = cl::sycl::buffer<buffer_data_type_t, 1>;
  auto thePointer = pMap.add_pointer(buffer_t(cl::sycl::range<1>{size}));
  
  return static_cast<void *>(thePointer);
}

template <bool ReUse = true, typename PointerMapper>
inline void SYCLfree(void *ptr, PointerMapper &pMap) {
  pMap.template remove_pointer<ReUse>(ptr);
}

template <typename PointerMapper>
inline void SYCLfreeAll(PointerMapper &pMap) {
  pMap.clear();
}

template <cl::sycl::access::mode AcMd, typename T>
struct RangeAccess {
  static const auto global_access = cl::sycl::access::target::global_buffer;
  static const auto is_place_holder = cl::sycl::access::placeholder::true_t;
  typedef T scalar_t;
  typedef scalar_t &ref_t;
  typedef typename cl::sycl::global_ptr<scalar_t>::pointer_t ptr_t;

  typedef cl::sycl::accessor<scalar_t, 1, AcMd, global_access, is_place_holder>
      accessor;

  typedef RangeAccess<AcMd, T> self_t;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE RangeAccess(accessor access,
                                                    size_t offset,
                                                    std::intptr_t virtual_ptr)
      : access_(access), offset_(offset), virtual_ptr_(virtual_ptr) {}

  RangeAccess(cl::sycl::buffer<scalar_t, 1> buff =
                  cl::sycl::buffer<scalar_t, 1>(cl::sycl::range<1>(1)))
      : access_{accessor{buff}}, offset_(0), virtual_ptr_(-1) {}

  RangeAccess(std::nullptr_t) : RangeAccess() {}
  
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ptr_t get_pointer() const {
    return (access_.get_pointer().get() + offset_);
  }
  template <typename Index>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE self_t &operator+=(Index offset) {
    offset_ += (offset);
    return *this;
  }
  template <typename Index>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE self_t operator+(Index offset) const {
    return self_t(access_, offset_ + offset, virtual_ptr_);
  }
  template <typename Index>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE self_t operator-(Index offset) const {
    return self_t(access_, offset_ - offset, virtual_ptr_);
  }
  template <typename Index>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE self_t &operator-=(Index offset) {
    offset_ -= offset;
    return *this;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE friend bool operator==(
      const RangeAccess &lhs, std::nullptr_t) {
    return ((lhs.virtual_ptr_ == -1));
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE friend bool operator!=(
      const RangeAccess &lhs, std::nullptr_t i) {
    return !(lhs == i);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE friend bool operator==(
      std::nullptr_t, const RangeAccess &rhs) {
    return ((rhs.virtual_ptr_ == -1));
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE friend bool operator!=(
      std::nullptr_t i, const RangeAccess &rhs) {
    return !(i == rhs);
  }
  
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE self_t &operator++() {
    offset_++;
    return (*this);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE self_t operator++(int i) {
    EIGEN_UNUSED_VARIABLE(i);
    self_t temp_iterator(*this);
    offset_++;
    return temp_iterator;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE std::ptrdiff_t get_size() const {
    return (access_.get_count() - offset_);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE std::ptrdiff_t get_offset() const {
    return offset_;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void set_offset(std::ptrdiff_t offset) {
    offset_ = offset;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ref_t operator*() const {
    return *get_pointer();
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ref_t operator*() {
    return *get_pointer();
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ptr_t operator->() = delete;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ref_t operator[](int x) {
    return *(get_pointer() + x);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ref_t operator[](int x) const {
    return *(get_pointer() + x);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE scalar_t *get_virtual_pointer() const {
    return reinterpret_cast<scalar_t *>(virtual_ptr_ +
                                        (offset_ * sizeof(scalar_t)));
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE explicit operator bool() const {
    return (virtual_ptr_ != -1);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE operator RangeAccess<AcMd, const T>() {
    return RangeAccess<AcMd, const T>(access_, offset_, virtual_ptr_);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  operator RangeAccess<AcMd, const T>() const {
    return RangeAccess<AcMd, const T>(access_, offset_, virtual_ptr_);
  }
  
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void bind(
      cl::sycl::handler &cgh) const {
    cgh.require(access_);
  }

 private:
  accessor access_;
  size_t offset_;
  std::intptr_t virtual_ptr_;  
};

template <cl::sycl::access::mode AcMd, typename T>
struct RangeAccess<AcMd, const T> : RangeAccess<AcMd, T> {
  typedef RangeAccess<AcMd, T> Base;
  using Base::Base;
};

}  
}  
}  

#endif  
