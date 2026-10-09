CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra

PREFIX  ?= /usr/local
BINDIR   = $(PREFIX)/bin

TARGET   = a119744

all: $(TARGET)

$(TARGET): a119744.c
	$(CC) $(CFLAGS) -o $@ a119744.c

install: $(TARGET)
	install -d $(BINDIR)
	install -m 755 $(TARGET) $(BINDIR)

clean:
	rm -f $(TARGET)

.PHONY: all install clean
