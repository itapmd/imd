//****************************************************************************
// Filename:            vvsocketio.h
// Project Affiliation: Virvo
// Funcionality:        Header file for vvsocketio.cpp
// Author:              Michael Poehnl
// Institution:         University of Stuttgart, HPC Center
// Operating Systems:   Linux, IRIX 6.5, Windows NT
// Creation Date:       01-10-31
//****************************************************************************

#ifndef _VVSOCKETIO_H_
#define _VVSOCKETIO_H_

#include "vvsocket.h"
#include "vvvoldesc.h"
#include "vvtoolshed.h"
#include "vvimage.h"
#include "vvvecmath.h"
#include "vvdebugmsg.h"

/** This class provides specific data transfer through sockets. 
  It requires the class vvSocket.<BR>
  Here is an example code fragment for a TCP sever which reads
  a volume from a socket and a TCP client which writes a volume to the
  socket:<BR>
  <PRE>
  
  Server:
  
  // Create a new TCP socket class instance, which is a server listening on port 17171
  
  vvSocketIO* sio = new vvSocketIO(17171 , vvSocket::VV_TCP);
  vvVolDesc* vd = new vvVolDesc();
  
  //Set the parameters of the socket( e.g. connection timer 3 sec, transfer
  timer 1.5 sec, socket buffer 65535 bytes, debuglevel 0)
  sio->set_sock_param(3.0f, 1.5f, 65535 , 0)
  
  // Get a volume 
  switch (sio->getVolume(vd)) 
  {
    case vvSocket::VV_OK:
      cerr << "Volume transferred successfully" << endl;
      break;
    case vvSocket::VV_ALLOC_ERROR:
      cerr << "Not enough memory" << endl;
      break;
    default:
      cerr << "Cannot read volume from socket" << endl;
      break;
  }
  delete sio;
  
  Client:
    
  // Create a new TCP socket class instance, which is a client and connects
  // to a server listening on port 17171
  
  char* servername = "buxdehude";
  vvSocketIO* sio = new vvSocketIO(17171 , servername, vvSocket::VV_TCP);
  vvVolDesc* vd = new vvVolDesc();

  //Set the parameters of the socket( e.g. connection timer 3 sec, transfer
  timer 1.5 sec, socket buffer 65535 bytes, debuglevel 0)
  sio->set_sock_param(3.0f, 1.5f, 65535 , 0);
  // Put a volume 
  switch (sio->putVolume(vd)) 
  {
    case vvSocket::VV_OK:
      cerr << "Volume transferred successfully" << endl;
      break;
    case vvSocket::VV_ALLOC_ERROR:
      cerr << "Not enough memory" << endl;
      break;
    default:
      cerr << "Cannot write volume to socket" << endl;
      break;
  }
  delete sio;  
  
  </PRE>
  @see vvSocket
  @author Michael Poehnl
*/
class vvSocketIO : public vvSocket
{
public:
        enum DataType     /// data type for get/putData
        {
                VV_UCHAR,
                VV_INT,
                VV_FLOAT
        };
        
        vvSocketIO(int, char*, vvSocket::SocketType, int clminport=0, int clmaxport=0);
        vvSocketIO(int, vvSocket::SocketType);
        ~vvSocketIO();
        ErrorType init();
        bool sock_action();
        ErrorType getVolume(vvVolDesc*);
        ErrorType putVolume(vvVolDesc*);  
        ErrorType getImage(vvImage*);
        ErrorType putImage(vvImage*);
        ErrorType getData(uchar**, int&);  //  unknown number and type
        ErrorType putData(uchar*, int);
        ErrorType getMatrix(vvMatrix4*);
        ErrorType putMatrix(vvMatrix4*);  
        ErrorType getData(void*, int, DataType); // known number and type
        ErrorType putData(void*, int, DataType);   
        void set_sock_param(float, float, int=65536, int=0);
};

#endif
