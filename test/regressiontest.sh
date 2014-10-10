#!/bin/bash

# Run from top-level match-vamp directory

VAMP_PATH=. \
    ../sonic-annotator/sonic-annotator \
    -d vamp:match-vamp-plugin:match \
    --multiplex \
    ~/Music/cc-kids-abrsm-dataset/ABRSM/Allegro\ in\ G.mp3 \
    ~/Music/cc-kids-abrsm-dataset/Kids/Allegro\ in\ G.mp3 \
    -w csv --csv-stdout 2>/dev/null | sed 's/^[^,]*,//' > /tmp/$$ || exit 1

diff /tmp/$$ `dirname $0`/expected.csv && echo Passed

rm /tmp/$$
