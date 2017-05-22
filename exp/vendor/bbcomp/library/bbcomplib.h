
#ifndef _bbcomplib_H_
#define _bbcomplib_H_


////////////////////////////////////////////////////////////
// DLL interface definition header file
//
// Include this header file when using the shared library.
// It provides the following macros:
//
// (*) DLL_DECLARATIONS
// (*) DLL_LOAD
// (*) DLL_LOAD_STATUS
// (*) DLL_CLOSE
//
// Use them as follows in your C/C++ program:
//
// #include <bbcomplib.h>
// DLL_DECLARATIONS
//
// int main()
// {
//   DLL_LOAD
//   if (DLL_LOAD_STATUS < 0) { printf("error loading library\n"); exit(EXIT_FAILURE); }
//   ...
//   DLL_CLOSE
// }
//
// After DLL_LOAD and before DLL_CLOSE, and as long as
// DLL_DECLARATIONS is in scope, the functions exported
// by the dynamic library can be called.
//


#include <assert.h>

// OS switches
#ifdef _WIN32
	#define _WINSOCKAPI_ 
	#include <windows.h>
	#define DLL_HANDLE_HELPER HINSTANCE libHandle;
	#define DLL_LOAD_HELPER libHandle = LoadLibrary("./bbcomp.dll"); if (! libHandle) { DLL_LOAD_STATUS = -1; }
	#define DLL_FUNCTION_HELPER(functiontype, name) name = (functiontype)(void*)GetProcAddress(libHandle, #name); if (! name) { DLL_LOAD_STATUS = -1; }
	#define DLL_CLOSE_HELPER FreeLibrary(libHandle);
#elif __linux__
	#include <dlfcn.h>
	#define DLL_HANDLE_HELPER void* libHandle;
	#define DLL_LOAD_HELPER libHandle = dlopen("./libbbcomp.so", RTLD_NOW); if (! libHandle) { DLL_LOAD_STATUS = -1; }
	#define DLL_FUNCTION_HELPER(functiontype, name) name = (functiontype)dlsym(libHandle, #name); if (! name) { DLL_LOAD_STATUS = -1; }
	#define DLL_CLOSE_HELPER dlclose(libHandle);
#elif __APPLE__
	#include <dlfcn.h>
	#define DLL_HANDLE_HELPER void* libHandle;
	#define DLL_LOAD_HELPER libHandle = dlopen("./libbbcomp.dylib", RTLD_NOW); if (! libHandle) { DLL_LOAD_STATUS = -1; }
	#define DLL_FUNCTION_HELPER(functiontype, name) name = (functiontype)dlsym(libHandle, #name); if (! name) { DLL_LOAD_STATUS = -1; }
	#define DLL_CLOSE_HELPER dlclose(libHandle);
#else
	#error unknown operating system
#endif


#ifdef __cplusplus
extern "C" {
#endif


// function pointer type definitions
typedef const char* stringtype;
typedef int (*int_func_int_string) (int, stringtype);
typedef int (*int_func) ();
typedef stringtype (*string_func_int) (int);
typedef int (*int_func_int) (int);
typedef int (*int_func_string) (stringtype);
typedef int (*int_func_string_string) (stringtype, stringtype);
typedef int (*int_func_point_point) (double*, double*);
typedef int (*int_func_int_point_point) (int, double*, double*);
typedef stringtype (*string_func) ();


#ifdef __cplusplus
} // extern "C"
#endif


#define DLL_DECLARATIONS \
	DLL_HANDLE_HELPER \
	int DLL_LOAD_STATUS; \
	int_func_int_string configure; \
	int_func_string_string login; \
	int_func numberOfTracks; \
	string_func_int trackName; \
	int_func_string setTrack; \
	int_func numberOfProblems; \
	int_func_int setProblem; \
	int_func dimension; \
	int_func numberOfObjectives; \
	int_func budget; \
	int_func evaluations; \
	int_func_point_point evaluate; \
	int_func_int_point_point history; \
	string_func errorMessage;

#define DLL_LOAD \
	DLL_LOAD_STATUS = 1; \
	DLL_LOAD_HELPER \
	DLL_FUNCTION_HELPER(int_func_int_string, configure) \
	DLL_FUNCTION_HELPER(int_func_string_string, login) \
	DLL_FUNCTION_HELPER(int_func, numberOfTracks) \
	DLL_FUNCTION_HELPER(string_func_int, trackName) \
	DLL_FUNCTION_HELPER(int_func_string, setTrack) \
	DLL_FUNCTION_HELPER(int_func, numberOfProblems) \
	DLL_FUNCTION_HELPER(int_func_int, setProblem) \
	DLL_FUNCTION_HELPER(int_func, dimension) \
	DLL_FUNCTION_HELPER(int_func, numberOfObjectives) \
	DLL_FUNCTION_HELPER(int_func, budget) \
	DLL_FUNCTION_HELPER(int_func, evaluations) \
	DLL_FUNCTION_HELPER(int_func_point_point, evaluate) \
	DLL_FUNCTION_HELPER(int_func_int_point_point, history) \
	DLL_FUNCTION_HELPER(string_func, errorMessage)

#define DLL_CLOSE \
	DLL_CLOSE_HELPER


#endif
