
phylip.o:     phylip.h

draw.o:   draw.c draw.h phylip.h 
	$(CC) $(DFLAGS) -c draw.c

draw2.o:   draw2.c draw.h phylip.h 
	$(CC) $(DFLAGS) -c draw2.c

drawgram.o:     drawgram.c draw.h phylip.h
	$(CC) $(DFLAGS) -c drawgram.c

drawgram:     drawgram.o draw.o draw2.o phylip.o
	$(CC) $(DFLAGS) draw.o draw2.o drawgram.o phylip.o -lm -o drawgram

all: drawgram

clean:
	@rm -f *.o