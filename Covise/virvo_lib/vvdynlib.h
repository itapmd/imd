//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef VV_DYNLIB_H
#define VV_DYNLIB_H

#ifdef WIN32
  #include <windows.h>
#endif

#ifdef __sgi
  #include <dlfcn.h>
#endif
#ifdef __hpux
  #include <dl.h>
#endif

#ifdef __hpux
  typedef shl_t VV_SHLIB_HANDLE;
#elif WIN32
  typedef HINSTANCE VV_SHLIB_HANDLE;
#else 
  typedef void *VV_SHLIB_HANDLE;
#endif
 
/** This class encapsulates the functionality of dynamic library loading.
  @author Uwe Woessner
*/
class vvDynLib
{
  public:
    static char* error(void);
    static VV_SHLIB_HANDLE open(const char* filename, int mode);
    static void* sym(VV_SHLIB_HANDLE handle, const char* symbolname);
    static int close(VV_SHLIB_HANDLE handle);
};

#endif
