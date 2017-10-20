CC = gcc
COPTS = -O2 # Use for distribution

INCLUDES = -Iinc

VERSION = -DVERSION=\"1.0.0\"
CFLAGS = $(COPTS) $(INCLUDES) $(VERSION)
LIBS =  -lm

OBJS =	basics.o hbind_itimer.o err_handle.o

libslideutils.a: $(OBJS)
	rm -f libhbindutils.a && \
        ar -r libhbindutils.a *.o && \
	ranlib libhbindutils.a

.c.o: 
	$(CC) $(CFLAGS) -c $< 

clean:
	rm -f *.a *.o *.gcda *.gcno 

touch:  
	touch *.c


