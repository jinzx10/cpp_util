#ifndef __MISCELLANEOUS_H__ 
#define __MISCELLANEOUS_H__

#include <cstdlib>
#include <sstream>

// read arguments from the command line
// no bound check!
template <int N = 1>
void readargs(char** args, std::string& var) {
	var = args[N];
}

template <int N = 1, typename T>
void readargs(char** args, T& var) {
	std::stringstream ss(args[N]);
	ss >> var;
}

template <int N = 1, typename T, typename ...Ts>
void readargs(char** args, T& var, Ts& ...rest) {
	readargs<N>(args, var);
	readargs<N+1>(args, rest...);
}

// mkdir
inline int mkdir(std::string const& dir) {
	std::string command = "mkdir -p " + dir;
	return std::system(command.c_str());
}

// touch
inline int touch(std::string const& file) {
	std::string command = "touch " + file;
	return std::system(command.c_str());
}

#endif
