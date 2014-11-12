#pragma once

/*
  This should be an almost drop in replacement for std::vector,
  with a few caveats.
  It allocates the first part of the space statically,
  and the rest in a std::vector, so if you need a bunch of these, and most
  will have bounded size, this should improve cache coherence substantially
  

*/
#include <iterator>
#include <type_traits>
#include <vector>
#include <cassert>

namespace benlib{

  template<typename T, size_t StaticSize = 8>
  class PreallocVector{

  public:
	
	//todo rename to "Iterator"
	template<typename VType, typename ContainerType>
	class PreallocVectorIterator : 
	  public std::iterator<std::random_access_iterator_tag,
						   VType>{
	
	public:
	  ~PreallocVectorIterator() = default; //makes clang whine
	  PreallocVectorIterator() = default; //don't do anything
	  PreallocVectorIterator(const PreallocVectorIterator&) = default;
	  PreallocVectorIterator& operator=(const PreallocVectorIterator&) = default;
	  //copy constructor and assignment

	  /*
	  template <typename otherVType, typename otherContainerType>
	  PreallocVectorIterator(const PreallocVectorIterator<otherVType, 
							 otherContainerType>& other)
		:index{other.index}, container{other.container}
	  {}
	  
	  template<typename otherVType, typename otherContainerType>
	  PreallocVectorIterator& 
	  operator=(const PreallocVectorIterator<otherVType, 
				otherContainerType>& other) noexcept{

		index = other.index;
		container = other.container;
		return *this;
	  }
	  */
	  //(in)equality comparison
	  template <typename otherVType, typename otherContainerType>
	  friend inline bool 
	  operator==(const PreallocVectorIterator<VType, ContainerType>& lhs,
				 const PreallocVectorIterator<otherVType, otherContainerType>& rhs) noexcept{
		
		return lhs.index == rhs.index;
	  }
	  template <typename otherVType, typename otherContainerType>
	  friend inline bool 
	  operator!=(const PreallocVectorIterator<VType, ContainerType>& lhs,
				 const PreallocVectorIterator<otherVType, otherContainerType>& rhs) noexcept{
		return !(lhs == rhs);
	  }
	  
	  template <typename otherVType, typename otherContainerType>
	  friend inline bool 
	  operator <(const PreallocVectorIterator<VType, ContainerType>& lhs,
				 const PreallocVectorIterator<otherVType, otherContainerType>& rhs) noexcept{
		
		return lhs.index < rhs.index;
	  }

	  template <typename otherVType, typename otherContainerType>
	  friend inline bool 
	  operator <=(const PreallocVectorIterator<VType, ContainerType>& lhs,
				  const PreallocVectorIterator<otherVType, otherContainerType>& rhs) noexcept{
		return lhs.index <= rhs.index;
	  }
	  
	  template <typename otherVType, typename otherContainerType>
	  friend inline bool 
	  operator >(const PreallocVectorIterator<VType, ContainerType>& lhs,
				 const PreallocVectorIterator<otherVType, otherContainerType>& rhs) noexcept{
		return lhs.index > rhs.index;
	  }
	  
	  template <typename otherVType, typename otherContainerType>
	  friend inline bool 
	  operator >=(const PreallocVectorIterator<VType, ContainerType>& lhs,
				  const PreallocVectorIterator<otherVType, otherContainerType>& rhs) noexcept{
		return lhs.index > rhs.index;
	  }



	  //increment
	  //template <typename intType>
	  PreallocVectorIterator& operator +=(std::ptrdiff_t amt) noexcept{
		//index += amt;
		index = static_cast<size_t>(static_cast<ptrdiff_t>(index) +
									amt);

		return *this;
	  }
	  PreallocVectorIterator& operator ++() noexcept{
		index++;
		return *this;
	  }
	  PreallocVectorIterator operator ++(int) noexcept{
		PreallocVectorIterator tmp{*this};
		operator++();
		return tmp;
	  }
	  
	  //decrement
	  //template <typename intType>
	  PreallocVectorIterator& operator -=(std::ptrdiff_t amt) noexcept{
		//index -= amount
		index = static_cast<size_t>(static_cast<ptrdiff_t>(index) -
									amt);
		return *this;
	  }
	  PreallocVectorIterator& operator --() noexcept{
		index--;
		return *this;
	  }
	  PreallocVectorIterator operator --(int) noexcept{
		PreallocVectorIterator tmp{*this};
		operator--();
		return tmp;
	  }
	  
	  //arithmetic
	  template<typename intType>
	  inline friend PreallocVectorIterator 
	  operator +(PreallocVectorIterator it, intType n){
		it += n;
		return it;
	  }
	  template<typename intType>
	  inline friend PreallocVectorIterator 
	  operator +(intType n, PreallocVectorIterator it){
		it += n;
		return it;
	  }
	  template<typename intType>
	  inline friend PreallocVectorIterator
	  operator -(PreallocVectorIterator it, intType n){
		it -= n;
		return it;
	  }

	  inline size_t
	  operator -(PreallocVectorIterator rhs)const {
		return index - rhs.index;
	  }


	  //dereferencing
	  const VType& operator*() const{
		return container->operator[](index);
	  }
	  const VType* operator->() const{
		return &(container->operator[](index));
	  }
	  const VType& operator[](size_t n) const{
		return container->operator[](index + n);
	  }

	  //non const derefs for nonconst iterator types
	  template<typename Dummy = VType>
	  typename std::enable_if<!std::is_const<Dummy>::value,
							  VType&>::type
	  operator*(){
		return container->operator[](index);
	  }
	  template<typename Dummy = VType>
	  typename std::enable_if<!std::is_const<Dummy>::value,
							  VType*>::type
	  operator->(){
		return &(container->operator[](index));
	  }
	  template<typename Dummy = VType>
	  typename std::enable_if<!std::is_const<Dummy>::value,
							  VType&>::type
	  operator[](size_t n){
		return container->operator[](index + n);
	  }
	  
	  
					 
	  template <typename U, size_t V>
	  friend class PreallocVector;
	private:
	  size_t index;
	  ContainerType* container;
	  PreallocVectorIterator(size_t _index, ContainerType* _container)
		:index{_index}, container{_container}
	  {}
	};

	
	

	using value_type = T;
	using reference = value_type&;
	using const_refernce = const reference;

	
	using iterator = PreallocVectorIterator<value_type, 
											PreallocVector<T, StaticSize>>;
	using const_iterator = 
	  PreallocVectorIterator<const value_type, 
							 const PreallocVector<T, StaticSize>>;
	


	//constructors
	PreallocVector()  //empty vector
	noexcept(std::is_nothrow_default_constructible<std::vector<T>>::value)
	  :size_{0}
	{}
	  
	// n default constructed objects
	explicit PreallocVector(size_t n)
	  :size_{n}
	{
	  for(size_t i = 0; i < n && i < StaticSize; ++i){
		new (data_() + i) value_type();
	  }
	  if(n > StaticSize){
		extraData_.resize(n - StaticSize);
	  }
	}
	
	// n not-default constructed objects
	PreallocVector(size_t n, const value_type& value)
	  : size_{n}
	{
	  for(size_t i = 0; i < n && i < StaticSize; ++i){
		new (data_() + i) value_type(value);
	  }
	  if(n > StaticSize){
		extraData_.assign(n - StaticSize, value);
	  }
	}
	
	// iterator range
	template <typename InIt>
	PreallocVector(InIt first, InIt last)
	  :size_{0}
	{
	  for(auto it = first; it != last; ++it){
		push_back(*it);
		++size_;
	  }
	}

	//copy
	template<size_t otherSize>
	PreallocVector(const PreallocVector<value_type, otherSize>& other)
	  :size_{0}
	{
	  assign(other.begin(), other.end());
	}

	

	inline void push_back(const value_type& val){
	  if(size_ < StaticSize){
		new (data_() + size_) value_type(val);
		++size_;
	  } else {
		extraData_.push_back(val);
		++size_;
	  }
	}
	
	inline void pop_back(){
	  assert(!empty());
	  if(size_ <= StaticSize){
		(data_() + size_ -1)->~value_type();
		--size_;
	  } else {
		extraData_.pop_back();
		--size_;
	  }
	}
	template<typename ... Args>
	inline void emplace_back(Args&&... args){
	  if(size_ < StaticSize){
		new(data_() + size_) value_type(std::forward<Args>(args)...);
		++size_;
	  } else {
		extraData_.emplace_back(std::forward<Args>(args)...);
		++size_;
	  }
	}

	inline void resize(size_t n){
	  if(n > size_){
		for(size_t i = size_; i < n; ++i){
		  emplace_back();
		}
	  } else {
		while(size_ > n){
		  pop_back();
		}
	  }
	}
	
	template <typename InIt>
	inline void assign(InIt first, InIt last){
	  size_t newSize = 0;
	  while(newSize < size_ && first != last){
		operator[](newSize) = *first++;
		++newSize;
	  }
	  while(first != last){
		push_back(*first++);
		++newSize;
	  }
	  while(newSize < size_){
		pop_back();
	  }
	  size_ = newSize;
	}

	inline void assign(size_t n, const value_type& val){
	  for(size_t i = 0; i < n && i < size_; ++i){
		operator[](i) = val;
	  }
	  for(size_t i = size_; i < n; ++i){
		push_back(val);
	  }
	  assert(size() == n);
	}


	inline value_type& operator[](size_t n){
	  if(n < StaticSize){
		return data_()[n];
	  }
	  return extraData_[n - StaticSize];
	}
	
	inline const value_type& operator[](size_t n) const {
	  if(n < StaticSize){
		return data_()[n];
	  }
	  return extraData_[n - StaticSize];
	}
	
	
	
	
	inline size_t size() const {return size_;}
	inline size_t staticSize() const { return StaticSize;}
	inline bool empty() const {return size_ == 0;}
	
	inline iterator begin(){ return iterator{0, this};}
	inline const_iterator begin() const {return const_iterator{0, this};}
	
	inline iterator end(){ return iterator{size_, this};}
	inline const_iterator end() const { return const_iterator{size_, this};}
		  
	inline value_type& front(){return data_()[0];}
	inline const value_type& front() const{return data_()[0];}
	
	inline value_type& back(){return operator[](size_ - 1);}
	inline const value_type& back() const {return operator[](size_ - 1);}

	//todo: insert operations


  private:
	size_t size_;
	char buffer_[sizeof(T)*StaticSize];
	//get at that buffer
	inline T* data_(){return reinterpret_cast<T*>(buffer_);} 
	inline const T* data_() const{return reinterpret_cast<const T*>(buffer_);} 

	std::vector<T> extraData_;


  };

}
