#!/usr/bin/env bash
# Run gmhmm-rs on every case in tests/cases.tsv and diff against tests/oracle.
# The date line in output headers is masked before comparing.
# Usage: tests/compare.sh [case-name-regex]
set -u
here=$(cd "$(dirname "$0")" && pwd)
# SUITE selects a case set: "." = gmhmme3 (default), "gmhmmp" = tests/gmhmmp
here=$here/${SUITE:-.}
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ "${SUITE:-.}" == "gmhmmp" ]]; then bin=${GMHMM_RS:-$root/target/release/gmhmmp-rs}; else bin=${GMHMM_RS:-$root/target/release/gmhmm-rs}; fi
ref=$here/oracle
tmp=$(mktemp -d)
pat=${1:-.}
pass=0; fail=0
mask() { sed -E 's/^(#+ ?|Date: |##date: )?(Mon|Tue|Wed|Thu|Fri|Sat|Sun) [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}$/<DATE>/' "$1"; }
cd "$here/data"
while IFS=$'\t' read -r name args; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    [[ "$name" =~ $pat ]] || continue
    # shellcheck disable=SC2086
    eval "set -- $args"
    "$bin" "$@" -o "$tmp/$name.out" > "$tmp/$name.stdout" 2> "$tmp/$name.stderr"
    echo $? > "$tmp/$name.rc"
    [[ -f "$name.b" ]] && mv "$name.b" "$tmp/$name.b"
    ok=1
    for ext in out stdout stderr rc b; do
        a=$ref/$name.$ext; b=$tmp/$name.$ext
        [[ -e "$a" || -e "$b" ]] || continue
        if ! diff -q <(mask "$a" 2>/dev/null) <(mask "$b" 2>/dev/null) > /dev/null; then
            ok=0
            echo "FAIL $name.$ext"
            diff <(mask "$a" 2>/dev/null) <(mask "$b" 2>/dev/null) | head -${DIFFLINES:-6}
        fi
    done
    if (( ok )); then pass=$((pass+1)); else fail=$((fail+1)); fi
done < "$here/cases.tsv"
echo "passed $pass, failed $fail"
rm -rf "$tmp"
(( fail == 0 ))
