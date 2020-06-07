#ifndef __STRING_HELPER_H__
#define __STRING_HELPER_H__

#include <cstdio>
#include <istream>
#include <iterator>
#include <string>
#include <cstring>
#include <type_traits>
#include <sstream>
#include <regex>
#include <iostream>
#include <complex>
#include "template_helper.h"

namespace cxut
{
	// remove certain leading and trailing characters
	inline std::string trim(std::string const& str, std::string rm = " \t") {
		auto start = str.find_first_not_of(rm);
		return (start == std::string::npos) ? 
			"" : str.substr(start, str.find_last_not_of(rm)-start+1);
	}
	
	// check if a string starts with a certain string
	inline bool start_with(std::string const& pre, std::string const& str, std::string skip = "") {
		auto start = str.find_first_not_of(skip);
		return (start == std::string::npos) ? (pre.size() ? false : true) : 
			(std::strncmp(pre.c_str(), str.substr(start).c_str(), pre.size()) == 0);
	}
	
	// replace the leading "~" of a string with ${HOME}
	inline std::string expand_leading_tilde(std::string const& dir, std::string skip = " \t") {
		auto start = dir.find_first_not_of(skip);
		return (start == std::string::npos || dir[start] != '~') ? dir :
			dir.substr(0, start) + std::getenv("HOME") + dir.substr(start+1);;
	}

	// remove all the occurences of certain characters from a string
	inline std::string remove_chars(std::string const& str, std::string const& rm = " ") {
		std::string result;
		std::copy_if(str.begin(), str.end(), std::back_inserter(result), [&rm] (char c) { return rm.find(c) == std::string::npos; });
		return result;
	}

	// regular expression
	// floating-point numbers
	const std::string fpn_regex = "([\\+-]?[0-9]*\\.?[0-9]*([eE][+-]?[0-9]+)?)";

	// c-style complex numbers: (real,imaginary) or (real)
	const std::string cx_fpn_regex = "(\\( *"+fpn_regex+" *, *"+fpn_regex+" *\\)|\\( *"+fpn_regex+" *\\))";

	// variant complex number format: a+bi/a+bj
	//const std::string cx_fpn_regex_var = "("+fpn_regex+" *(\\+|-) *"+fpn_regex+"(i|j))";

	// all supported numbers
	const std::string num_regex = "(("+fpn_regex+")|("+cx_fpn_regex+"))";
	//const std::string num_regex = "(("+fpn_regex+")|("+cx_fpn_regex+")|("+cx_fpn_regex_var+"))";

	// brace-enclosed list
	inline std::string bracelist_regex(std::string const& entry) {
		return "((\\{ *\\})|(\\{ *"+entry+"( *, *"+entry+")* *\\}))";
	}


	inline int conv_string(std::string const& from, bool& to) {
		std::stringstream ss(from);
		ss >> to;
		if (ss.fail()) {
			ss.clear();
			ss >> std::boolalpha >> to;
		}
		return ss.good() ? 0 : 1;
	}

	template <typename T>
	typename std::enable_if< std::is_trivial<T>::value || is_std_complex<T>::value, int>::type conv_string(std::string const& from, T& to) {
		std::stringstream ss(from);
		ss >> to;
		return ss.good() ? 0 : 1;
	}

	inline int conv_string(std::string const& from, std::string& to) {
		to = from;
		return 0;
	}


	inline int is_balanced(std::string const& str, std::string const& balance = "{}()[]") {
		std::stack<char> stk;
		if (balance.size()%2) {
			std::cerr << "is_balanced: error: characters to balance must be paired." << std::endl;
			exit(EXIT_FAILURE);
		}

		std::string left, right;
		for (size_t i = 0; i != balance.size()/2; ++i) {
			left += balance[2*i];
			right += balance[2*i+1];
		}

		for (size_t i = 0; i != str.size(); ++i) {
			if ( left.find(str[i]) != std::string::npos ) {
				stk.push(str[i]);
				continue;
			}
			
			auto idx = right.find(str[i]);
			if ( idx != std::string::npos ) {
				if ( stk.empty() || stk.top() != left[idx] )
					return false;
				else
					stk.pop();
			}
		}
		return stk.empty();
	}

	
	inline std::istream& getuntil(std::istream& ss, std::string& str, std::string const& delim="\n", std::string const& balance="{}()[]") {
		str.clear();
		if (balance.size()%2) {
			std::cerr << "getuntil: error: characters to balance must be paired." << std::endl;
			exit(EXIT_FAILURE);
		}

		std::stack<char> stk;
		bool is_balanced = true;
		std::string left, right;
		for (size_t i = 0; i != balance.size()/2; ++i) {
			left += balance[2*i];
			right += balance[2*i+1];
		}

		auto reach_delim = [&ss, &delim] () -> bool { 
			return delim.find(ss.peek()) != std::string::npos; 
		};

		auto reach_space = [&ss] () -> bool { 
			return ss.peek() == ' '; 
		};

		auto reach_end = [&ss] () -> bool { 
			return ss.peek() == std::char_traits<char>::eof(); 
		};

		auto ignore_delim_and_space = [&] () {
			while ( !reach_end() && ( reach_delim() || reach_space() )  )
				ss.ignore(1);
		};

		// ignore leading delimiters and whitespaces
		ignore_delim_and_space();

		while ( !reach_end() ) {
			//std::cout << "next char = " << ss.peek() << std::endl;
			if ( left.find(ss.peek()) != std::string::npos ) {
				stk.push(ss.peek());
				is_balanced = false;
			} else {
				auto idx = right.find(ss.peek());
				if ( idx != std::string::npos ) {
					if ( stk.empty() || stk.top() != left[idx] ) {
						std::cerr << "getuntil: error: input is not balanced." << std::endl;
						exit(EXIT_FAILURE);
					} else {
						stk.pop();
						is_balanced = stk.empty();
					}
				}
			}

			if ( !reach_delim() || !is_balanced ) {
				str += ss.get();
			} else {
				ignore_delim_and_space();
				break;
			}
		}
		return ss;
	}

}

#endif
