#ifndef _CLIENT_H_
#define _CLIENT_H_

int initClient(int *soc, int i, char  *server_name);
void connect_server(char token);
double get_maximum(double *array, int size);
double get_minimum(double *array, int size);
void scale_array(double *array, int size, double min, double max, double *buffer);
void write_bytes_to_file(char filename[256],char *byte_array, int size);

#endif
