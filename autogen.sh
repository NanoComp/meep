#!/bin/sh

echo "arg: $1"

touch hsrc/.depend

if test "x$1" = "x--verbose"; then
    verbose=yes
else
    verbose=no
fi

if test $verbose = yes; then
    autoreconf --verbose --install --symlink
else
    autoreconf --install --symlink
fi

config=good # hackery so darcs_test still outputs config.log w/failed configure

./configure --enable-maintainer-mode || config=bad

if test $verbose = yes; then
    cat config.log
fi

test $config = bad && exit 1

(cd hsrc; make depend)

