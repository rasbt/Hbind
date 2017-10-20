#include <sys/time.h>
#include <signal.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "hbind_itimer.h"

static long real_elapsed_sec;
static long virt_elapsed_sec;
static long prof_elapsed_sec;
static long timer_interval_sec;


void hbind_itimer_sig_handler(int signum, siginfo_t *siginfo, void *ucontext)
{
  if(signum == SIGALRM) real_elapsed_sec += timer_interval_sec;
  else if(signum == SIGVTALRM) virt_elapsed_sec += timer_interval_sec;
  else if(signum == SIGPROF) prof_elapsed_sec += timer_interval_sec;
}


int hbind_itimer_start(int interval_sec_in)
{
  struct sigaction alrm_act;
  struct itimerval tval;
  timer_interval_sec = interval_sec_in;
  alrm_act.sa_sigaction = hbind_itimer_sig_handler;
  alrm_act.sa_flags = SA_SIGINFO | SA_RESTART;
  sigemptyset(&(alrm_act.sa_mask));
  real_elapsed_sec = 0;
  virt_elapsed_sec = 0;
  prof_elapsed_sec = 0;
  tval.it_value.tv_sec = timer_interval_sec;
  tval.it_value.tv_usec = 0;
  tval.it_interval.tv_sec = timer_interval_sec;
  tval.it_interval.tv_usec = 0;

  if(sigaction(SIGALRM, &alrm_act, NULL)){
    fprintf(stderr, "Unable to create SIGALRM for the profile timer");
    return -1;
  }else if(sigaction(SIGVTALRM, &alrm_act, NULL)){
    fprintf(stderr, "Unable to create SIGVTALRM for the profile timer");
    return -1;
  }else if(sigaction(SIGPROF, &alrm_act, NULL)){
    fprintf(stderr, "Unable to create SIGPROF for the profile timer");
    return -1;
  }else if( setitimer(ITIMER_REAL, &tval, 0) == -1){
    fprintf(stderr, "Unable to set the real timer: %s", strerror(errno)); 
    return -1;
  }else if( setitimer(ITIMER_VIRTUAL, &tval, 0) == -1){
    fprintf(stderr, "Unable to set the virtual timer: %s", strerror(errno)); 
    return -1;
  }else if( setitimer(ITIMER_PROF, &tval, 0) == -1){
    fprintf(stderr, "Unable to set the profile timer: %s", strerror(errno)); 
    return -1;
  }
  return 0;
}

int get_itimer(double *rv, int which, long elapsed_sec)
{
  struct itimerval tval;
  if( getitimer(which, &tval) == -1){
    fprintf(stderr, "Unable to get the current value of the timer: %s",
            strerror(errno));
    *rv = -1.0;
    return -1;
  }
  /* timers count down from starting value*/
  *rv = elapsed_sec + timer_interval_sec - tval.it_value.tv_sec;
  *rv -= tval.it_value.tv_usec / 1000000.0;
  if(*rv < 0) *rv = 0;
  return 0;

}

int hbind_itimer_get(double *real, double *virt, double *prof)
{
  if (get_itimer(real, ITIMER_REAL, real_elapsed_sec) == 0
      && get_itimer(virt, ITIMER_VIRTUAL, virt_elapsed_sec) == 0  
      &&  get_itimer(prof, ITIMER_PROF, prof_elapsed_sec) == 0 ) 
    return 0;
  else
    return 1;
}

char* hbind_get_local_time(char *buf, size_t buflen)
{
  struct timeval tv;
  struct timezone tz;
  struct tm ltime;

  if(buflen <= 0) return buf;
  if(gettimeofday(&tv, &tz) == -1){
    fprintf(stderr, "Unable to get time using the \"gettimeofday\" function: %s",
            strerror(errno));
    buf[0] = 0; 
    return buf;
  }

  if(!localtime_r(&(tv.tv_sec), &ltime)){
    fprintf(stderr, "Unable to get the local time using \"localtime_r\": %s",
            strerror(errno));
    buf[0] = 0; 
    return buf;
  }

  if(my_strftime(buf, buflen, "%c", &ltime) == 0) buf[buflen - 1] = 0;
  return buf;
}

size_t my_strftime(char *s, size_t max, const char  *fmt,
                   const struct tm *tm)
{ return strftime(s, max, fmt, tm); }
