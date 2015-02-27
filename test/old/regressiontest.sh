#!/bin/bash

mydir=$(dirname "$0")

sonic-annotator --minversion 1.1 || exit 1

VAMP_PATH="${mydir}/../.." \
    sonic-annotator \
    -d vamp:match-vamp-plugin:match \
    --multiplex \
    ~/Music/cc-kids-abrsm-dataset/ABRSM/Allegro\ in\ G.mp3 \
    ~/Music/cc-kids-abrsm-dataset/Kids/Allegro\ in\ G.mp3 \
    -w csv --csv-stdout 2>/dev/null | sed 's/^[^,]*,//' > /tmp/$$ || exit 1

sdiff -w 78 /tmp/$$ `dirname $0`/expected.csv && echo Passed

rm /tmp/$$
