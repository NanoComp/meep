#!/bin/sh

test -r /usr/share/misc/config.sub && \
   cp -f /usr/share/misc/config.sub config.sub
test -r /usr/share/misc/config.guess && \
   cp -f /usr/share/misc/config.guess config.guess
autoheader
autoconf
