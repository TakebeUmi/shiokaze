/*
**	main.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Feb 25, 2018.
**
**	Permission is hereby granted, free of charge, to any person obtaining a copy of
**	this software and associated documentation files (the "Software"), to deal in
**	the Software without restriction, including without limitation the rights to use,
**	copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
**	Software, and to permit persons to whom the Software is furnished to do so,
**	subject to the following conditions:
**
**	The above copyright notice and this permission notice shall be included in all copies
**	or substantial portions of the Software.
**
**	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
**	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
**	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
**	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
**	CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
**	OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//
#include <iostream>
#include <dlfcn.h>
#include <thread>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include <shiokaze/core/filesystem.h>
//
SHKZ_USING_NAMESPACE
//
// https://stackoverflow.com/questions/77005/how-to-automatically-generate-a-stacktrace-when-my-program-crashes
static void handler (int sig) {
	//
	void *array[10];
	size_t size;
	//
	// get void*'s for all entries on the stack
	size = backtrace(array, 10);
	//
	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	exit(1);
}
//
int main ( int argc, const char* argv[] ) {
	//
	signal(SIGSEGV, handler); // install our handler
	pthread_mutex_t mutex;
	pthread_mutex_init(&mutex, NULL);
	//
	std::thread dummy([](){}); dummy.join(); // Dummy code to enforce linking against thread
	const auto handle = ::dlopen(filesystem::resolve_libname("shiokaze_ui").c_str(),RTLD_LAZY);
	int result (0);
	if( ! handle ) {
		std::cout << "Could not open the library: " << ::dlerror() << std::endl;
	} else {
		const auto run_func = ::dlsym(handle,"run");
		if( ! run_func ) {
			std::cout << "Could not load the function: " << ::dlerror() << std::endl;
		} else {
			const auto func = reinterpret_cast<int(*)(int argc, const char* argv[])>(run_func);
			result = func(argc,argv);
			if(::dlclose(handle)) {
				std::cout << "Could not close the handle: " << ::dlerror() << std::endl;
			}
		}
	}
	//
	return result;
}
//