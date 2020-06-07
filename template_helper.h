#ifndef __TEMPLATE_HELPER_H__
#define __TEMPLATE_HELPER_H__

#include <type_traits>
#include <complex>

namespace cxut {

	template <typename ...> using void_t = void;

	template <typename T, typename = void>
	struct is_std_complex : std::false_type {};

	template <typename T>
	struct is_std_complex<T, typename std::enable_if<std::is_same<std::complex<typename T::value_type>, T>::value, void>::type> : std::true_type {};


	template <typename F, typename ...Args>
	class is_invocable
	{
	    typedef char yes;
	    typedef char no[2];
	
		template <typename F1, typename ...Args1>
		using return_t = decltype( std::declval<F1>()( std::declval<Args1>()... ) ); // remove reference?
	
	    template <typename F1, typename ...Args1>
	    static yes& test(return_t<F1, Args1...>* );
	
	    template <typename F1, typename ...Args1>
	    static no& test(...);
	
		template <bool, typename, typename ...>
		struct try_return { using type = void; };
		
		template <typename F1, typename ...Args1>
		struct try_return<true, F1, Args1...> { using type = return_t<F1, Args...>; };
	
	
		public:
	
	    static const bool value = ( sizeof(test<F, Args...>(nullptr)) == sizeof(yes) );
		using return_type = typename try_return<value, F, Args...>::type;
	
	};

}


#endif
