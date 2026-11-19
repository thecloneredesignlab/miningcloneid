/* Read isolation for the future-evaluator interface test, not scientific code. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <unistd.h>
static int deny(const char *path){
 if(!path)return 0;
 char full[PATH_MAX],root[PATH_MAX];const char *r=getenv("CG_PRIVATE_ROOT");if(!r)return 0;
 if(!realpath(path,full))return 0;
 snprintf(root,sizeof(root),"%s/inputs/private",r);
 int blocked=!strncmp(full,root,strlen(root)) || strstr(full,"/20261116_culture_step05_bridge_design_01/")!=NULL;
 const char *tasks[]={"I01","I01_R1","P01","R01","R02","Q01","Q02"};
 for(int i=0;i<7;i++){snprintf(root,sizeof(root),"%s/verification/%s/",r,tasks[i]);if(!strncmp(full,root,strlen(root)))blocked=1;}
 /* I01's source-verified guard library itself is loaded before this check can run. */
 if(blocked){const char *log=getenv("CG_READ_DENIAL_LOG");if(log){int fd=syscall(SYS_openat,AT_FDCWD,log,O_WRONLY|O_CREAT|O_APPEND,0600);if(fd>=0){syscall(SYS_write,fd,full,strlen(full));syscall(SYS_write,fd,"\n",1);syscall(SYS_close,fd);}}errno=EACCES;return 1;}
 return 0;
}
#define OPEN(N) int N(const char*p,int f,...){va_list a;va_start(a,f);mode_t m=(f&O_CREAT)?va_arg(a,int):0;va_end(a);if(deny(p))return -1;int(*fn)(const char*,int,...)=dlsym(RTLD_NEXT,#N);return fn(p,f,m);}
OPEN(open)
OPEN(open64)
#define FOPEN(N) FILE*N(const char*p,const char*m){if(deny(p))return NULL;FILE*(*fn)(const char*,const char*)=dlsym(RTLD_NEXT,#N);return fn(p,m);}
FOPEN(fopen)
FOPEN(fopen64)
