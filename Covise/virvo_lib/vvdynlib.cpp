//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifdef _STANDARD_C_PLUS_PLUS
  #include <iostream>
  using std::cerr;
  using std::endl;
#else
  #include <iostream.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __linux
  #include <string.h>
  #include <dlfcn.h>
#endif

#include "vvdynlib.h"

int vvDynLib::close(VV_SHLIB_HANDLE handle)
{
#ifdef WIN32
  FreeLibrary (handle);
  return 1;
  
#elif __hpux
  // HP-UX 10.x and 32-bit 11.00 do not pay attention to the ref count when
  // unloading a dynamic lib.  So, if the ref count is more than 1, do not
  // unload the lib.  This will cause a library loaded more than once to
  // not be unloaded until the process runs down, but that's life.  It's
  // better than unloading a library that's in use.
  // So far as I know, there's no way to decrement the refcnt that the kernel
  // is looking at - the shl_descriptor is a copy of what the kernel has, not
  // the actual struct.
  // On 64-bit HP-UX using dlopen, this problem has been fixed.
  struct shl_descriptor  desc;
  if (shl_gethandle_r(handle, &desc) == -1)
    return -1;
  if (desc.ref_count > 1)
    return 1;
#if defined(__GNUC__) || __cplusplus >= 199707L
  shl_unload(handle);
#else
  cxxshl_unload(handle);
#endif  /* aC++ vs. Hp C++ */
  return 1;
  
#else

#ifdef __sun4
  // SunOS4 does not automatically call _fini()!
  void *ptr;
  ptr = sym(handle, "_fini");

  if (ptr != 0)
    (*((int (*)(void)) ptr)) (); // Call _fini hook explicitly.
#endif

  dlclose(handle);
  return 1;
#endif 
}

char* vvDynLib::error()
{
#ifdef __hpux
  return strerror(errno);
#elif WIN32
  static char buf[128];
  FormatMessageA (FORMAT_MESSAGE_FROM_SYSTEM,
                  NULL,
                  GetLastError (),
                  0,
                  buf,
                  sizeof buf,
                  NULL);
  return buf;
#else
  return (char*)dlerror();
#endif 
}

VV_SHLIB_HANDLE vvDynLib::open(const char* filename, int mode)
{
  void* handle;
  
#ifdef WIN32
  handle = LoadLibraryA (filename);
#elif __hpux
#if defined(__GNUC__) || __cplusplus >= 199707L
  handle = shl_load(filename, mode, 0L);
#else
  handle = cxxshl_load(filename, mode, 0L);
#endif  
#else
  handle = dlopen(filename, mode);
#endif

  if(handle == NULL)
  {
    cerr << error() << endl;
  }

#ifdef __sun4
  if (handle != 0)
    {
      void *ptr;
      // Some systems (e.g., SunOS4) do not automatically call _init(), so
      // we'll have to call it manually.

      ptr = sym(handle, "_init");

      if (ptr != 0 && (*((int (*)(void)) ptr)) () == -1) // Call _init hook explicitly.
        {
          // Close down the handle to prevent leaks.
          close(handle);
          return 0;
        }
    }
#endif

#ifdef WIN32
  mode = mode;    // prevent warning
#endif

  return (VV_SHLIB_HANDLE)handle;
}

void* vvDynLib::sym(VV_SHLIB_HANDLE handle, const char* symbolname)
{
#ifdef WIN32
  return GetProcAddress(handle, symbolname);
#elif __hpux
  void *value;
  int status;
  shl_t _handle = handle;
  status = shl_findsym(&_handle, symbolname, TYPE_UNDEFINED, &value);
  return status == 0 ? value : NULL;
#else
  return dlsym (handle, symbolname);
#endif
}
