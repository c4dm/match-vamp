#!/bin/bash

# Test that aligning A-B gives the same results as aligning B-A

mydir=$(dirname "$0")

sonic-annotator --minversion 1.1 || exit 1

file1="$1"
file2="$2"

usage() {
    echo "Usage: $0 a.wav b.wav"
    exit 2
}

test -n "$file1" || usage
test -f "$file1" || usage
test -n "$file2" || usage
test -f "$file2" || usage

VAMP_PATH="${mydir}/../.." \
    sonic-annotator \
    -d vamp:match-vamp-plugin:match:a_b \
    --multiplex \
    "$file1" \
    "$file2" \
    -w csv --csv-stdout 2>/dev/null | sed 's/^"[^"]*"//' > /tmp/1_$$ || exit 1

VAMP_PATH="${mydir}/../.." \
    sonic-annotator \
    -d vamp:match-vamp-plugin:match:b_a \
    --multiplex \
    "$file2" \
    "$file1" \
    -w csv --csv-stdout 2>/dev/null | sed 's/^"[^"]*"//' > /tmp/2_$$ || exit 1

sdiff -w 78 /tmp/1_$$ /tmp/2_$$ && echo Passed

rm /tmp/1_$$ /tmp/2_$$
