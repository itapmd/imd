#ifndef _SOCKUTIL_H_
#define _SOCKUTIL_H_
#include <sys/types.h>
int WriteFull(int filedes, const void *buffer, int nbytes);
int ReadFull(int filedes, const void *buffer, int nbytes);
unsigned long GetIP(const char *name);
int OpenServerSocket(u_short MyPort);
int OpenNBServerSocket(u_short MyPort);
int OpenClientSocket(u_long toIP, u_short toPort);
void WriteSync(int fd);
void ReadSync(int fd);
#endif
