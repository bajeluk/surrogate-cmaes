
#ifndef _bbcomplib_H_
#define _bbcomplib_H_


#include <assert.h>


////////////////////////////////////////////////////////////
// OS switches for DLL init / quit / export
//
#ifdef _WIN32
	
	#define DLL_INIT static void initializer(void)
	#define DLL_QUIT static void finalizer(void)
#elif __linux__
	#define DLL_EXPORT __attribute__((visibility("default")))
	#define DLL_INIT __attribute__((constructor)) static void initializer(void)
	#define DLL_QUIT __attribute__((destructor)) static void finalizer(void)
#elif __APPLE__
	#define DLL_EXPORT __attribute__((visibility("default")))
	#define DLL_INIT __attribute__((constructor)) static void initializer(void)
	#define DLL_QUIT __attribute__((destructor)) static void finalizer(void)
#else
	#error unknown operating system
#endif


// OS switches
#ifdef _WIN32
	#define _WINSOCKAPI_ 
	#include <windows.h>
	#define DLL_HANDLE_HELPER HINSTANCE libHandle;
	#define DLL_LOAD_HELPER libHandle = LoadLibrary("bbcomp.dll"); assert(libHandle);
	#define DLL_FUNCTION_HELPER(functiontype, name) name = (functiontype)(void*)GetProcAddress(libHandle, #name); assert(name);
	#define DLL_CLOSE_HELPER FreeLibrary(libHandle);
#elif __linux__
	#include <dlfcn.h>
	#define DLL_HANDLE_HELPER void* libHandle;
	#define DLL_LOAD_HELPER libHandle = dlopen("libbbcomp.dylib", RTLD_NOW); assert(libHandle);
	#define DLL_FUNCTION_HELPER(functiontype, name) name = (functiontype)dlsym(libHandle, #name); assert(name);
	#define DLL_CLOSE_HELPER dlclose(libHandle);
#elif __APPLE__
	#include <dlfcn.h>
	#define DLL_HANDLE_HELPER void* libHandle;
	#define DLL_LOAD_HELPER libHandle = dlopen("libbbcomp.dylib", RTLD_NOW); assert(libHandle);
	#define DLL_FUNCTION_HELPER(functiontype, name) name = (functiontype)dlsym(libHandle, #name); assert(name);
	#define DLL_CLOSE_HELPER dlclose(libHandle);
#else
	#error unknown operating system
#endif


#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
	#define DLL_EXPORT __declspec( dllexport ) 
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


DLL_EXPORT int configure(int history, const char* logfilepath);
DLL_EXPORT int login(stringtype username, stringtype password);
DLL_EXPORT int numberOfTracks();
DLL_EXPORT stringtype trackName(int trackindex);
DLL_EXPORT int setTrack(stringtype trackname);
DLL_EXPORT int numberOfProblems();
DLL_EXPORT int setProblem(int problemID);
DLL_EXPORT int dimension();
DLL_EXPORT int numberOfObjectives();
DLL_EXPORT int budget();
DLL_EXPORT int evaluations();
DLL_EXPORT int evaluate(double* point, double* value);
DLL_EXPORT int history(int index, double* point, double* value);
DLL_EXPORT stringtype errorMessage();

#ifdef __cplusplus
} // extern "C"
#endif


#define DLL_DECLARATIONS \
	DLL_HANDLE_HELPER \
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
