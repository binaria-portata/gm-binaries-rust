#!/usr/bin/env bash
# Generate reference outputs with the original gmhmme3 (GeneMark-ETP build).
# Usage: GMHMME3_REF=/path/to/gmhmme3 tests/make_oracle.sh
# Writes tests/oracle/<case>.{out,stdout,stderr,rc}. The case list is shared
# with tests/compare.sh through tests/cases.tsv.
set -u
here=$(cd "$(dirname "$0")" && pwd)
# SUITE selects a case set: "." = gmhmme3 (default), "gmhmmp" = tests/gmhmmp
here=$here/${SUITE:-.}
ref=${GMHMME3_REF:?set GMHMME3_REF to the original gmhmme3 binary}
out=$here/oracle
mkdir -p "$out"
cd "$here/data"
while IFS=$'\t' read -r name args; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    # shellcheck disable=SC2086
    eval "set -- $args"
    "$ref" "$@" -o "$out/$name.out" > "$out/$name.stdout" 2> "$out/$name.stderr"
    echo $? > "$out/$name.rc"
    # side files written by -b go to <case>.b
    [[ -f "$name.b" ]] && mv "$name.b" "$out/$name.b"
    echo "$name rc=$(cat "$out/$name.rc")"
done < "$here/cases.tsv"
