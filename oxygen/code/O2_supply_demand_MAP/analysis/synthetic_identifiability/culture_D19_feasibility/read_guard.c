/* Scoped read barrier for observation-free workers. */
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
static int prefix(const char*x,const char*y){return !strncmp(x,y,strlen(y));}
static int deny(const char*path){
 if(!path)return 0;const char*r=getenv("REC_STUDY_ROOT");if(!r)return 0;
 char full[PATH_MAX],root[PATH_MAX],buf[PATH_MAX];if(!realpath(path,full)||!realpath(r,root))return 0;
 const char*mode=getenv("REC_MODE");const char*origin=getenv("REC_ORIGIN");int diag=mode&&!strcmp(mode,"DIAG");
 snprintf(buf,sizeof(buf),"%s/verification/P0/feasibility_contract.rds",root);if(!strcmp(full,buf))return 0;
 snprintf(buf,sizeof(buf),"%s/verification/P0/read_guard.so",root);if(!strcmp(full,buf))return 0;
 if(origin){snprintf(buf,sizeof(buf),"%s/verification/P0/rng/%s.rds",root,origin);if(!strcmp(full,buf))return 0;}
 if(diag){snprintf(buf,sizeof(buf),"%s/verification/P0/old_units.rds",root);if(!strcmp(full,buf))return 0;}
 int blocked=0;
 snprintf(buf,sizeof(buf),"%s/inputs/private/",root);if(prefix(full,buf))blocked=1;
 snprintf(buf,sizeof(buf),"%s/inputs/diagnostic/",root);if(prefix(full,buf)&&!diag)blocked=1;
 snprintf(buf,sizeof(buf),"%s/verification/P0/",root);if(prefix(full,buf))blocked=1;
 snprintf(buf,sizeof(buf),"%s/verification/P1/",root);if(prefix(full,buf))blocked=1;
 snprintf(buf,sizeof(buf),"%s/verification/P1_",root);if(prefix(full,buf))blocked=1;
 snprintf(buf,sizeof(buf),"%s/reports/",root);if(prefix(full,buf))blocked=1;
 int own=!strcmp(full,root)||(prefix(full,root)&&full[strlen(root)]=='/');
 if(strstr(full,"/oxygen/results/")&&!own)blocked=1;
 snprintf(buf,sizeof(buf),"%s/verification/F_",root);if(prefix(full,buf)){
  if(!origin)blocked=1;else{snprintf(buf,sizeof(buf),"%s/verification/%s_",root,origin);if(!prefix(full,buf))blocked=1;}
 }
 if(!blocked)return 0;const char*log=getenv("REC_READ_DENIAL_LOG");if(log){int fd=syscall(SYS_openat,AT_FDCWD,log,O_WRONLY|O_CREAT|O_APPEND,0600);if(fd>=0){syscall(SYS_write,fd,full,strlen(full));syscall(SYS_write,fd,"\n",1);syscall(SYS_close,fd);}}
 errno=EACCES;return 1;
}
static mode_t modearg(int f,va_list a){return (f&O_CREAT)||((f&O_TMPFILE)==O_TMPFILE)?(mode_t)va_arg(a,int):0;}
#define OPEN(N) int N(const char*p,int f,...){va_list a;va_start(a,f);mode_t m=modearg(f,a);va_end(a);if(deny(p))return -1;int(*fn)(const char*,int,...)=dlsym(RTLD_NEXT,#N);return fn(p,f,m);}
OPEN(open)
OPEN(open64)
#define FOPEN(N) FILE*N(const char*p,const char*m){if(deny(p))return NULL;FILE*(*fn)(const char*,const char*)=dlsym(RTLD_NEXT,#N);return fn(p,m);}
FOPEN(fopen)
FOPEN(fopen64)
#define OPENAT(N) int N(int d,const char*p,int f,...){va_list a;va_start(a,f);mode_t m=modearg(f,a);va_end(a);char q[PATH_MAX];const char*x=p;if(p&&p[0]!='/'&&d!=AT_FDCWD){char proc[100];snprintf(proc,sizeof(proc),"/proc/self/fd/%d",d);ssize_t k=readlink(proc,q,sizeof(q)-1);if(k>=0){q[k]=0;size_t z=strlen(q);snprintf(q+z,sizeof(q)-z,"/%s",p);x=q;}}if(deny(x))return -1;int(*fn)(int,const char*,int,...)=dlsym(RTLD_NEXT,#N);return fn(d,p,f,m);}
OPENAT(openat)
OPENAT(openat64)
