#ifndef __PARSER_H__
#define __PARSER_H__

/* A keyword parser 
 *
 * Usage: 
 * Parser p({"key1", "key2", ...}); 
 * p.parse(filename); // filename must be a std::string
 * p.pour(val1, val2, ...); // key1-val1, key2-val2, etc.   */


#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "string_helper.h"


namespace cxut {

	struct Parser
	{
		Parser(std::vector<std::string> const& keys_) : keys(keys_), vals(keys.size()) {}
	
		std::vector<std::string> keys;
		std::vector<std::string> vals;
	
		void reset(std::vector<std::string> const& keys_ = {});
		void parse(std::string const& file);

		template <typename T> int assign(std::string const& key, T& val);
		template <size_t N = 0, typename T> void pour(T&);
		template <size_t N = 0, typename T, typename ...Rs> void pour(T&, Rs&...);
	
		private:
		void size_check(size_t);
	};


	inline void Parser::reset(std::vector<std::string> const& keys_) {
		keys = keys_;
		vals = std::vector<std::string>(keys.size(), "");
	}

	inline void Parser::parse(std::string const& file) {
		std::fstream fs(file);
		std::string str;
		while( getuntil(fs, str, ",;\n", "{}()[]") ) {
			str = trim(str);
			for (size_t i = 0; i != keys.size(); ++i) {
				if (start_with(keys[i], str)) {
					str.erase(0, keys[i].size());
					vals[i] = trim(str);
					break;
				}
			}
		}	
	}

	template <typename T>
	int Parser::assign(const std::string& key, T& val) {
		for (size_t i = 0; i != keys.size(); ++i) {
			if (key.compare(keys[i]) == 0) {
				conv_string(vals[i], val);
				return 0;
			}
		}
		std::cerr << "key \"" << key << "\" not found." << std::endl;
		return 1;
	}

	template <size_t N, typename T>
	void Parser::pour(T& val) {
		size_check(N);
		conv_string(vals[N], val);
	}
	
	template <size_t N, typename T, typename ...Rs>
	void Parser::pour(T& val, Rs& ...args) {
		pour<N>(val);
		pour<N+1>(args...);
	}

	inline void Parser::size_check(size_t N) {
		if ( N >= keys.size() ) {
			std::cerr << "Parser error: too many variables: expect " << keys.size() << " or less." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

}


#endif
