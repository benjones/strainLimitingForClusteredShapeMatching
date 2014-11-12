#pragma once

//#include <type_traits>

namespace benlib{



  //inspired by cppitertools
  template <typename T>
  //for integer types with step size 1, we can be really fast
  class FastRange{
  public:
	class Iterator :
	  public std::iterator<std::random_access_iterator_tag, T> {
	  
	private:
	  T value;
	public:
	  //default constructors, operators, etc
	  friend inline bool operator ==(const Iterator& lhs, const Iterator& rhs){
		return lhs.value == rhs.value;
	  }
	  friend inline bool operator !=(const Iterator& lhs, const Iterator& rhs){
		return lhs.value != rhs.value;
	  }
	  friend inline bool operator <(const Iterator& lhs, const Iterator& rhs){
		return lhs.value < rhs.value;
	  }
	  friend inline bool operator <=(const Iterator& lhs, const Iterator& rhs){
		return lhs.value <= rhs.value;
	  }
	  friend inline bool operator >(const Iterator& lhs, const Iterator& rhs){
		return lhs.value > rhs.value;
	  }
	  friend inline bool operator >=(const Iterator& lhs, const Iterator& rhs){
		return lhs.value >= rhs.value;
	  }
	  Iterator& operator ++() noexcept{
		++value;
		return *this;
	  }
	  Iterator operator ++(int) noexcept{
		Iterator tmp{*this};
		operator ++();
		return tmp;
	  }
	  Iterator& operator --() noexcept{
		--value;
		return *this;
	  }
	  Iterator operator --(int) noexcept{
		Iterator tmp{*this};
		operator --();
		return tmp;
	  }
	  template <typename IntType>
	  Iterator& operator +=(IntType n) noexcept{
		value += n;
		return *this;
	  }
	  template <typename IntType>
	  Iterator& operator -=(IntType n) noexcept{
		value -= n;
		return *this;
	  }
	  template<typename IntType>
	  inline friend Iterator operator +(Iterator lhs, IntType n) noexcept{
		lhs += n;
		return lhs;
	  }
	  template<typename IntType>
	  inline friend Iterator operator -(Iterator lhs, IntType n) noexcept{
		lhs -= n;
		return lhs;
	  }
	  inline friend ptrdiff_t
	  operator -(Iterator lhs, Iterator rhs) noexcept{
		return lhs.value - rhs.value;
	  }

	  inline const T& operator *() const {
		return value;
	  }
	  inline T* operator ->() const {
		return &value;
	  }

	  T operator[](size_t n) const{
		Iterator tmp{*this};
		tmp += n;
		return *tmp;
	  }
	  //no real advantage for making this private...
	  Iterator(T _value) noexcept
	  : value{_value}
	  {}
	};

	using iterator = Iterator;
	using const_iterator = Iterator;

	FastRange(T _begin, T _end) noexcept
	:begin_{_begin}, end_{_end}
	{}

	
	FastRange(T _end) noexcept
	:FastRange{0, _end}
	{ }


	
	Iterator begin() const {return begin_;}
	Iterator end() const {return end_;}
  private:
	const Iterator begin_, end_;

  };

  template<typename T>
  class Range{
  public:
	
	//Range will only return const iterators(doesn't make sense to write to an iterator)
	class Iterator : 
	  public std::iterator<std::random_access_iterator_tag, T> {
	private:
	  T value;
	  T step; //you could want an unsigned iterator with a negative step...
	  bool isEnd;

	  bool comesBefore(const Iterator& other, std::true_type /*T is unsigned */) const{
		return value < other.value;
	  }
	  bool comesBefore(const Iterator& other, std::false_type /*T is signed */) const{
		return (step < 0) ? (value > other.value) : (value < other.value);
	  }
	  
	  
	public:
	  //default constructors/destructors
	  
	  friend inline bool operator==(const Iterator& lhs, const Iterator& rhs){
		return (lhs.isEnd && rhs.isEnd) ||
		  lhs.value == rhs.value;
	  }
	  //this is tricky because we don't want to overshoot end
	  // and stl uses != even when it should use < 
	  friend inline bool operator!=(const Iterator& lhs, const Iterator& rhs){
		return (lhs.isEnd && rhs.isEnd) ? false : //ends are ==
		  (rhs.isEnd ? lhs.comesBefore(rhs, typename std::is_unsigned<T>::type()) : //if one is end
		   (lhs.isEnd ? rhs.comesBefore(lhs, typename std::is_unsigned<T>::type()) : //true if "in bounds"
			lhs.value != rhs.value)); //if neither are ends, just do !=

	  }
	  //these comparisons don't care about whether the iterator is the end or not
	  friend inline bool operator <(const Iterator& lhs, const Iterator& rhs){
		assert((lhs.step < 0 && rhs.step < 0) ||
			   (lhs.step > 0 && rhs.step > 0));
		return lhs.step > 0 ? (lhs.value < rhs.value) : (lhs.value > rhs.value);
	  }
	  friend inline bool operator <=(const Iterator& lhs, const Iterator& rhs){
		assert((lhs.step < 0 && rhs.step < 0) ||
			   (lhs.step > 0 && rhs.step > 0));
		return lhs.step > 0 ? (lhs.value <= rhs.value) : (lhs.value >= rhs.value);
	  }
	  friend inline bool operator <(const Iterator& lhs, const Iterator& rhs){
		assert((lhs.step < 0 && rhs.step < 0) ||
			   (lhs.step > 0 && rhs.step > 0));
		return lhs.step > 0 ? (lhs.value > rhs.value) : (lhs.value < rhs.value);
	  }
	  friend inline bool operator <(const Iterator& lhs, const Iterator& rhs){
		assert((lhs.step < 0 && rhs.step < 0) ||
			   (lhs.step > 0 && rhs.step > 0));
		return lhs.step > 0 ? (lhs.value >= rhs.value) : (lhs.value <= rhs.value);
	  }
	  
	  Iterator& operator++() noexcept{
		value += step;
		return *this;
	  }
	  Iterator operator++(int) noexcept{
		Iterator tmp{*this};
		operator++();
		return tmp;
	  }
	  
	  Iterator& operator--() noexcept{
		value -= step;
		return *this;
	  }

	  Iterator operator--(int) noexcept{
		Iterator tmp{*this};
		operator--();
		return tmp;
	  }

	  template<typename IntType>
	  Iterator& operator +=(IntType n) noexcept{
		value += step*n;
		return *this;
	  }
	  
	  template<typename IntType>
	  Iterator& operator -=(IntType n) noexcept{
		value -= step*n;
		return *this;
	  }
	  
	  template<typename IntType>
	  inline friend Iterator
	  operator +(Iterator lhs, IntType n){
		lhs += n;
		return lhs;
	  }
	  template<typename IntType>
	  inline friend Iterator
	  operator -(Iterator lhs, IntType n){
		lhs -= n;
		return lhs;
	  }

	  inline ptrdiff_t 
	  operator -(const Iterator& rhs) const{
		assert(rhs.step == step); //distance doesn't make sense otherwise
		return static_cast<ptrdiff_t>((value - rhs.value)/step);
	  }

	  const T& operator*() const{
		return value;
	  }
	  const T* operator->() const{
		return &value; //is this useful?
	  }
	  T operator[](size_t n) const{
		Iterator tmp{*this};
		tmp += n;
		return tmp.value;
	  }
	  
	  friend class Range;
	private:
	  Iterator(T _value, T _step, bool _isEnd) noexcept
	  :value{_value}, step{_step}, isEnd{_isEnd}
	  {}
	};

	using iterator = Iterator;
	using const_iterator = Iterator;

	//using the python conventions for constructing these
	Range(T _begin, T _end, T step) noexcept
	:begin_{_begin, step, false},
	  end_{_end, step, true}
	{}
	Range(T _begin, T _end) noexcept
	: Range{_begin, _end, T{1}}
	{}
	Range(T _end) noexcept
	:Range{0, _end, 1}
	{}
	
	
	Iterator begin() const {return begin_;}
	Iterator end() const {return end_;}
  private:
	const Iterator begin_, end_;

  };

  //again, using the python conventions
  template<typename T>
  Range<T> range(T start, T stop, T step){ return Range<T>{start, stop, step};}

  template<typename T>
  typename std::conditional<std::is_integral<T>::value,
							FastRange<T>,
							Range<T> >::type
  range(T start, T stop){return typename std::conditional<std::is_integral<T>::value,
														  FastRange<T>,
														  Range<T> >::type{start, stop};}
  
  template<typename T>
  typename std::conditional<std::is_integral<T>::value,
							FastRange<T>,
							Range<T> >::type
  range(T stop){return typename std::conditional<std::is_integral<T>::value,
												 FastRange<T>,
												 Range<T> >::type{stop};}
  
}
