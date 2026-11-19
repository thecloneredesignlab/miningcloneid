/* Canonical role allowlists; no private response or foreign-origin handoff. */
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
#include <sys/stat.h>
#include <unistd.h>
static int under(const char*x,const char*y){size_t n=strlen(y);return !strncmp(x,y,n)&&(x[n]=='/'||x[n]==0);}
static int deny(const char*path){
 const char*r=getenv("REC_STUDY_ROOT"),*m=getenv("REC_MODE"),*o=getenv("REC_ORIGIN");
 if(!path||!r||!m||!o)return 0;
 char full[PATH_MAX],root[PATH_MAX],buf[PATH_MAX];if(!realpath(path,full)||!realpath(r,root))return 0;
 struct stat st;if(stat(full,&st)||!S_ISREG(st.st_mode))return 0;
 int own=under(full,root),allowed=0,blocked=0;
 if(own){
  const char *suffix=full+strlen(root);
  if(!strncmp(suffix,"/release/",9)||!strncmp(suffix,"/runtime/",9)||!strncmp(suffix,"/inputs/production/",19))allowed=1;
  if(!strcmp(suffix,"/inputs/public/contract.rds")||!strcmp(suffix,"/inputs/public/parameter_roles.tsv"))allowed=1;
  if(!strcmp(suffix,"/verification/T00/read_guard.so"))allowed=1;
  snprintf(buf,sizeof(buf),"/verification/T00/rng/%s.rds",o);if(!strcmp(suffix,buf))allowed=1;
  snprintf(buf,sizeof(buf),"/verification/T00/identities/%s.rds",o);if(!strcmp(suffix,buf))allowed=1;
  snprintf(buf,sizeof(buf),"/verification/%s_",o);if(!strncmp(suffix,buf,strlen(buf)))allowed=1;
  /* Canary output belongs to its own allocation. */
  const char *out=getenv("CULTURE_ALLOWED_WRITE_ROOT");if(out&&under(full,out))allowed=1;
  if(!strcmp(m,"D")){
   if(!strcmp(suffix,"/inputs/public/D19.rds"))allowed=1;
   snprintf(buf,sizeof(buf),"/inputs/banks/%s.rds",o);if(!strcmp(suffix,buf))allowed=1;
   snprintf(buf,sizeof(buf),"/controller/F_%c_bank_path.txt",o[2]);if(!strcmp(suffix,buf))allowed=1;
   /* Validated reserve bank remains in the same family, separate from other D origins. */
   if(strlen(o)==3){snprintf(buf,sizeof(buf),"/verification/F_%c_",o[2]);if(!strncmp(suffix,buf,strlen(buf)))allowed=1;}
  }
  if(!strcmp(m,"F")){
   snprintf(buf,sizeof(buf),"/controller/%s_bank_path.txt",o);if(!strcmp(suffix,buf))allowed=1;
  }
  blocked=!allowed;
 }else{
  if(strstr(full,"/oxygen/results/"))blocked=1;
  const char*cut=strstr(root,"/oxygen/results/");if(cut){char project[PATH_MAX];size_t n=cut-root;memcpy(project,root,n);project[n]=0;if(under(full,project))blocked=1;}
 }
 const char*b=strrchr(full,'/');b=b?b+1:full;
 if(!strcmp(b,"fit_data.Rds")||!strcmp(b,"jobs_2N.Rds")||!strcmp(b,"jobs_4N.Rds"))blocked=1;
 if(!blocked)return 0;
 const char*log=getenv("REC_READ_DENIAL_LOG");if(log){int fd=syscall(SYS_openat,AT_FDCWD,log,O_WRONLY|O_CREAT|O_APPEND,0600);if(fd>=0){syscall(SYS_write,fd,full,strlen(full));syscall(SYS_write,fd,"\n",1);syscall(SYS_close,fd);}}
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
