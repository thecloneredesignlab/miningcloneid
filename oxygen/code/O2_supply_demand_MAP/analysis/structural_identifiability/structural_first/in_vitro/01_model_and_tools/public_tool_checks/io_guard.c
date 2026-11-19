/* Runtime policy guard only. Does not implement or modify a scientific method. */
#define _GNU_SOURCE
#include <ctype.h>
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <unistd.h>

static int forbidden(const char *path) {
    if (!path) return 0;
    char lower[8192]; size_t i;
    for(i=0; path[i] && i<sizeof(lower)-1; ++i) lower[i]=(char)tolower((unsigned char)path[i]);
    lower[i]=0;
    const char *base=strrchr(lower,'/'); base=base?base+1:lower;
    if(strstr(lower,".json") || strstr(lower,"libjvm") || strstr(lower,"libjava.") ||
       strstr(lower,"libjson") || strstr(lower,"libmozjs") || strstr(lower,"libv8")) return 1;
    size_t n=strlen(base);
    return !strcmp(base,"java") || !strcmp(base,"node") || !strcmp(base,"nodejs") ||
           !strcmp(base,"js") || (n>=3 && !strcmp(base+n-3,".js"));
}

static int reject(const char *kind,const char *path) {
    if(!forbidden(path)) return 0;
    const char *dest=getenv("CULTURE_GUARD_LOG");
    if(dest) {
        int fd=(int)syscall(SYS_openat,AT_FDCWD,dest,O_WRONLY|O_CREAT|O_APPEND,0600);
        if(fd>=0) {
            char line[9000]; int n=snprintf(line,sizeof(line),"BLOCKED\t%s\t%s\n",kind,path);
            if(n>0) syscall(SYS_write,fd,line,(size_t)n<sizeof(line)?(size_t)n:sizeof(line)-1);
            syscall(SYS_close,fd);
        }
    }
    errno=EACCES; return 1;
}

static mode_t getmode(int flags,va_list ap) {
    return (flags&O_CREAT) || ((flags&O_TMPFILE)==O_TMPFILE) ? (mode_t)va_arg(ap,int) : 0;
}

#define OPEN_WRAPPER(NAME) \
int NAME(const char *path,int flags,...) { \
    va_list ap; va_start(ap,flags); mode_t mode=getmode(flags,ap); va_end(ap); \
    if(reject(#NAME,path)) return -1; \
    int (*real)(const char*,int,...)=dlsym(RTLD_NEXT,#NAME); return real(path,flags,mode); }
OPEN_WRAPPER(open)
OPEN_WRAPPER(open64)

#define OPENAT_WRAPPER(NAME) \
int NAME(int dirfd,const char *path,int flags,...) { \
    va_list ap; va_start(ap,flags); mode_t mode=getmode(flags,ap); va_end(ap); \
    if(reject(#NAME,path)) return -1; \
    int (*real)(int,const char*,int,...)=dlsym(RTLD_NEXT,#NAME); return real(dirfd,path,flags,mode); }
OPENAT_WRAPPER(openat)
OPENAT_WRAPPER(openat64)

#define FOPEN_WRAPPER(NAME) \
FILE *NAME(const char *path,const char *mode) { \
    if(reject(#NAME,path)) return NULL; \
    FILE *(*real)(const char*,const char*)=dlsym(RTLD_NEXT,#NAME); return real(path,mode); }
FOPEN_WRAPPER(fopen)
FOPEN_WRAPPER(fopen64)

#define FREOPEN_WRAPPER(NAME) \
FILE *NAME(const char *path,const char *mode,FILE *stream) { \
    if(reject(#NAME,path)) return NULL; \
    FILE *(*real)(const char*,const char*,FILE*)=dlsym(RTLD_NEXT,#NAME); return real(path,mode,stream); }
FREOPEN_WRAPPER(freopen)
FREOPEN_WRAPPER(freopen64)

int creat(const char *path,mode_t mode) {
    if(reject("creat",path)) return -1;
    int (*real)(const char*,mode_t)=dlsym(RTLD_NEXT,"creat"); return real(path,mode);
}
void *dlopen(const char *path,int flags) {
    if(reject("dlopen",path)) return NULL;
    void *(*real)(const char*,int)=dlsym(RTLD_NEXT,"dlopen"); return real(path,flags);
}
int execve(const char *path,char *const argv[],char *const env[]) {
    if(reject("execve",path)) return -1;
    int (*real)(const char*,char *const[],char *const[])=dlsym(RTLD_NEXT,"execve"); return real(path,argv,env);
}
