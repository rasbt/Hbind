
#ifndef _TIMER_H_INCLUDED
#define _TIMER_H_INCLUDED
#include <time.h>

int hbind_itimer_start(int interval_sec_in);

int hbind_itimer_get(double *real, double *virt, double *prof);

char* hbind_get_local_time(char *buf, size_t buflen);

/*! Wrapper to eliminate buggy gcc warning about using "%c" in strftime.*/
/*! From the man page for strftime:
 * Some buggy versions of gcc complain about the use of %c: warning: `%c'
 * yields only last 2 digits of year in some locales.  Of course
 * programmers are encouraged to use %c, it gives the preferred date
 * and time representation. One meets all kinds of strange obfuscations to
 * circumvent this gcc problem. A relatively clean one is to add this
 * intermediate function.
 *
 * @param s Character array of size max
 * @param max Size of the array s
 * @param fmt Format string
 * @param tm The tm structure to convert to a string.
 * @return The number of chars written to s if the time fits, else 0.
 */
size_t my_strftime(char *s, size_t max, const char  *fmt, const struct tm *tm);

#endif
