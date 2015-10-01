/* -------------------------------------------------
   All binary files are in C format
   Fortran does not have a binary format by default
   
   The following functions handle
   - open
   - write
   - read
   - close
   ------------------------------------------------ */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include <sys/stat.h>
#include <sys/types.h>

/* Constants */
/*const static int size_file = 64;*/
const static int none = 0;
const static int little = 1;
const static int big = 2;

/* Global variables */
static FILE* fp[10];
static int used[10] = {0,0,0,0,0,0,0,0,0,0};

/* Definition of the functions*/
void BINARY_FILE_OPEN(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2);
void BINARY_FILE_CLOSE(int* unit, int* ierr);
void BINARY_FILE_WRITE(int* unit, const void* buffer, int* count, int* size, int* ierr);
void BINARY_FILE_READ(int* unit, void* buffer, int* count, int* size, int* ierr);
void CREATE_FOLDER(const char* file, const int size);

/* Fortran Wrapper*/
void binary_file_open(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2){
  BINARY_FILE_OPEN(unit, file, mode, ierr, size1, size2);
}
void binary_file_open_(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2){
  BINARY_FILE_OPEN(unit, file, mode, ierr, size1, size2);
}
void binary_file_open__(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2){
  BINARY_FILE_OPEN(unit, file, mode, ierr, size1, size2);
}
void binary_file_close(int* unit, int* ierr){
  BINARY_FILE_CLOSE(unit, ierr);
}
void binary_file_close_(int* unit, int* ierr){
  BINARY_FILE_CLOSE(unit, ierr);
}
void binary_file_close__(int* unit, int* ierr){
  BINARY_FILE_CLOSE(unit, ierr);
}
void binary_file_write(int* unit, const void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_file_write_(int* unit, const void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_file_write__(int* unit, const void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_WRITE(unit, buffer, count, size, ierr);
}
void binary_file_read(int* unit, void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_READ(unit, buffer, count, size, ierr);
}
void binary_file_read_(int* unit, void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_READ(unit, buffer, count, size, ierr);
}
void binary_file_read__(int* unit, void* buffer, int* count, int* size, int* ierr){
  BINARY_FILE_READ(unit, buffer, count, size, ierr);
}
void create_folder(const char* file, const int size){
  CREATE_FOLDER(file,size);
}
void create_folder_(const char* file, const int size){
  CREATE_FOLDER(file,size);
}
void create_folder__(const char* file, const int size){
  CREATE_FOLDER(file,size);
}


/* Open the file*/
void BINARY_FILE_OPEN(int* unit, const char* file, const char* mode, int* ierr, const int size1, const int size2){

  char *tmp;
  int test;

  for( *unit=0; *unit<10 && used[*unit]==1; (*unit)++) {}
  if (*unit==10) {
    printf("BINARY_FILE_OPEN : Too many files open\n");
    *ierr = 1;
    return;
  }
  used[*unit] = 1;
    
  char filename[64/*size_file*/];
  strncpy(filename,file,(size_t) size1);
  filename[size1] = '\0';

  char format[64/*size_file*/];
  strncpy(format,mode,(size_t) size2);
  format[size2] = '\0';

  fp[*unit] = fopen(filename,format);

  if (fp[*unit]==NULL){
    printf("BINARY_FILE_OPEN : Error openning the file\n");    
    *ierr = 2;
    return;
  }else
    *ierr = 0;

  return;
}

/* Close the file*/
void BINARY_FILE_CLOSE(int* unit, int* ierr){

  fclose(fp[*unit]);
  used[*unit] = 0;
  return;
}

/* Write to a file*/
void BINARY_FILE_WRITE(int* unit, const void* buffer, int* count, int* size,  int* ierr){

  int i,n;
  unsigned char* tmp;
  char * pos;

  tmp = malloc(*size);
  pos = (char*) buffer;

  for (n=0; n<*count; n++){
    memcpy(tmp,pos,*size);
    for (i=0; i<*size; i++){
      fputc(tmp[i], fp[*unit]);
    }
    pos += *size;
  }
  
  free(tmp);
  *ierr = 0;
  return;
}

/* Read from a file*/
void BINARY_FILE_READ(int* unit, void* buffer, int* count, int* size,  int* ierr){

  int i,n;
  unsigned char* tmp;
  char * pos;

  tmp = malloc(*size);
  pos = (char*) buffer;
  
  for (n=0; n<*count; n++){
    for (i=0; i<*size; i++){
      tmp[i] = fgetc(fp[*unit]);
    }
    memcpy(pos,tmp,*size);
    pos += *size;
  }
  
  free(tmp);
  *ierr = 0;
  return;
}

/* Create a folder */
void CREATE_FOLDER(const char* file, const int size){
  
  char dirname[64];
  strncpy(dirname,file,(size_t) size);
  dirname[size] = '\0';
  
  int status;
  status = mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}
