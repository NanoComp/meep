#!/bin/sh

touch hsrc/.depend

autoreconf --verbose --install --symlink

./configure --enable-maintainer-mode

(cd hsrc; make depend)

