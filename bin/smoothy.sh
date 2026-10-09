#!/bin/sh
set -e

scriptdir=$(dirname "$0")

# Starts the smoothy module: all of Genestrip's and ft's goals plus the smoothy ones. Unlike
# bin/genestrip.sh, it runs lib/smoothy.jar, which is built by `mvn -pl smoothy -am package`.
# Heap as in bin/genestrip.sh; override per run with GS_XMX if a machine has more or less.
: "${GS_XMX:=56g}"

java -Xmx"$GS_XMX" -jar $scriptdir/../lib/smoothy.jar -d $scriptdir/../data "$@"
