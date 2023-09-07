CC = gcc
CFLAGS = -Wall -Werror -O3

SRCS = bshuf_h5filter.c bitshuffle.c bitshuffle_core.c iochain.c
OBJS = $(SRCS:.c=.o)
TARGET = libbitshuffle.a

all: $(TARGET)

$(TARGET): $(OBJS)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
