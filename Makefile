objects: heat_serial

CC = g++
DEBUG = -g
CFLAGS = -Wall $(DEBUG)


heat_serial : heat_serial.cpp
	$(CC) $(CFLAGS) heat_serial.cpp -o $(@)

clean:
	\rm *.o heat_serial



