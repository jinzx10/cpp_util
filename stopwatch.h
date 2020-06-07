#ifndef __STOPWATCH_H__
#define __STOPWATCH_H__

#include <chrono>
#include <unordered_map>
#include <iostream>
#include <string>
#include "template_helper.h"
#include <type_traits>

namespace cxut {

	class Stopwatch
	{
	public:

		Stopwatch(): t_start(), dur_store(), is_ticking() {}
	
		void start(size_t const& i = 0);
		void pause(size_t const& i = 0);
		void report(size_t const& i = 0);
		void reset(size_t const& i = 0);

		template <typename F, typename ...Args>
		typename std::enable_if<cxut::is_invocable<F, Args...>::value, void>::type timeit(size_t const& N, F f, Args const& ...args);

		template <typename F, typename ...Args>
		typename std::enable_if<cxut::is_invocable<F, Args...>::value, void>::type timeit(F f, Args const& ...args) { timeit(10ul, f, args...); };
	
		void timeit(...) { std::cerr << "Stopwatch error: timeit: invalid call." << std::endl; }


	private:

		using iclock		= std::chrono::high_resolution_clock;
		using dur_t			= std::chrono::duration<double>;
	
		std::unordered_map<size_t, iclock::time_point> t_start;
		std::unordered_map<size_t, dur_t> dur_store;
		std::unordered_map<size_t, bool> is_ticking;
	
		int get_status(size_t const& i) { return (is_ticking.find(i) == is_ticking.end()) ? -1 : is_ticking[i]; }

		template <typename F, typename ...Args>
		int try_it(size_t const& N, F f, Args const& ...args) { return N ? f(args...), try_it(N-1, f, args...) : 0; }
	};


	inline void Stopwatch::start(size_t const& i) {
		switch (get_status(i)) {
			case 1: std::cerr << "Stopwatch #" << i << " is already ticking." << std::endl;
					break;
			case -1: dur_store[i] = dur_t::zero(); // fallthrough
			case 0: t_start[i] = iclock::now(); 
					is_ticking[i] = true;
		}
	}

	inline void Stopwatch::pause(size_t const& i) { 
		switch (get_status(i)) {
			case -1: std::cerr << "Stopwatch #" << i << " does not exist." << std::endl;
					 break;
			case 0: std::cerr << "Stopwatch #" << i << " is already paused." << std::endl;
					break;
			case 1: dur_store[i] += iclock::now() - t_start[i]; 
					is_ticking[i] = false;
		}
	}

	inline void Stopwatch::report(size_t const& i) { 
		switch (get_status(i)) {
			case -1: std::cerr << "Stopwatch #" << i << " does not exist." << std::endl;
					 break;
			case 0: std::cout << "elapsed time = " << dur_store[i].count() << " seconds (stopwatch #" << i << ", paused)." << std::endl; 
					break;
			case 1: dur_t dur = dur_store[i] + static_cast<dur_t>(iclock::now() - t_start[i]); 
					std::cout << "elapsed time = " << dur.count() << " seconds (stopwatch #" << i << ", ticking)." << std::endl; 
		}
	}

	inline void Stopwatch::reset(size_t const& i) { 
		dur_store[i] = dur_store[i].zero();
		is_ticking[i] = false;
	}

	template <typename F, typename ...Args>
	typename std::enable_if<is_invocable<F, Args...>::value, void>::type Stopwatch::timeit(size_t const& N, F f, Args const& ...args) {
		dur_t dur;
		iclock::time_point start = iclock::now();
		try_it<F, Args...>(N, f, args...);
		dur = iclock::now() - start;
		std::cout << "average elapsed time for " << N << " trials = " << dur.count() / ( N ? N : 1 ) << " seconds." << std::endl;
	}
}



#endif
