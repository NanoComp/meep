#!/bin/sh
#
# Build a version number, possibly including the Git revision number. Called
# from configure.ac via AC_INIT.
#

#set -x

VERSION=$1
TAG=$2

# Include the Git revision when the release tag is "alpha"
if [ "$TAG" = "alpha" -a -e .git ]; then
    REV=$(git rev-parse --short=6 HEAD)
fi

# Print the version string, depending on how many components are available
if [ "$TAG" != "" -a "$REV" != "" ]; then
    printf "%s-%s.%s" $VERSION $TAG $REV
elif [ "$TAG" != "" ]; then
    printf "%s-%s" $VERSION $TAG
else
    printf "%s" $VERSION
fi
