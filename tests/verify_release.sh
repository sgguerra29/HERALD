#!/bin/bash
# Verify that this release reproduces output from a previous HERALD run.
#
# Both stages are compared against files you already have, so this checks the
# whole pipeline rather than any single function.
#
# IMPORTANT: pass the same --mmf and -s you used originally. The default
# changed to 1 - 1/s in this release, so omitting --mmf will legitimately
# produce different results and the comparison will be meaningless.
#
# Usage:
#   ./verify_release.sh -s 10 --mmf 0.55 \
#       --sam        original_alignment.sam \
#       --fragments  previous_fragments.fasta \
#       --realigned  realigned_fragments.sam \
#       --results    previous_results.txt
#
# --fragments and --results are each optional; whichever you supply is compared.

set -u

HERALD_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SIZE=10
MMF=""
SAM=""
PREV_FRAGS=""
REALIGNED=""
PREV_RESULTS=""
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--size)       SIZE="$2"; shift 2 ;;
    --mmf)           MMF="$2"; shift 2 ;;
    --sam)           SAM="$2"; shift 2 ;;
    --fragments)     PREV_FRAGS="$2"; shift 2 ;;
    --realigned)     REALIGNED="$2"; shift 2 ;;
    --results)       PREV_RESULTS="$2"; shift 2 ;;
    *)               EXTRA+=("$1"); shift ;;
  esac
done

if [[ -z "$MMF" ]]; then
  echo "Refusing to run without --mmf."
  echo "The default changed in this release; pass the value your original run used."
  exit 2
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp -r "$HERALD_ROOT/src" "$WORK/src"   # run in a temp tree, leave your files alone

FAIL=0

hr() { printf '%s\n' "----------------------------------------------------------"; }

if [[ -n "$SAM" ]]; then
  hr; echo "STAGE 1  fragmenting $SAM"; hr
  python3 "$WORK/src/launcher.py" -i "$SAM" -s "$SIZE" --mmf "$MMF" "${EXTRA[@]+"${EXTRA[@]}"}" || exit 1
  NEW_FRAGS="$WORK/tmp/fragments.fasta"
  echo "  fragments written: $(grep -c '^>' "$NEW_FRAGS")"

  if [[ -n "$PREV_FRAGS" ]]; then
    if diff -q <(sort "$NEW_FRAGS") <(sort "$PREV_FRAGS") >/dev/null; then
      echo "  IDENTICAL to $PREV_FRAGS"
    else
      echo "  DIFFERS from $PREV_FRAGS"
      echo "    previous: $(grep -c '^>' "$PREV_FRAGS") fragments"
      echo "    first differing records:"
      diff <(grep '^>' "$NEW_FRAGS" | sort) <(grep '^>' "$PREV_FRAGS" | sort) | head -10
      FAIL=1
    fi
  fi
fi

if [[ -n "$REALIGNED" ]]; then
  hr; echo "STAGE 2  scoring $REALIGNED"; hr
  python3 "$WORK/src/launcher.py" -i "$REALIGNED" -s "$SIZE" --mmf "$MMF" \
      --results -o verify "${EXTRA[@]+"${EXTRA[@]}"}" || exit 1
  NEW_RESULTS="$WORK/Output/Post_Proc_Fragment_Results_verify.txt"
  [[ -f "$NEW_RESULTS" ]] || NEW_RESULTS="$WORK/Output/Fragment_Results_verify.txt"
  echo "  candidate reads: $(grep -c '_frag_' "$NEW_RESULTS" || true) fragment rows"

  if [[ -n "$PREV_RESULTS" ]]; then
    if diff -q "$NEW_RESULTS" "$PREV_RESULTS" >/dev/null; then
      echo "  IDENTICAL to $PREV_RESULTS"
    else
      echo "  DIFFERS from $PREV_RESULTS"
      echo "    candidate reads, new vs previous:"
      echo "      $(grep -c '_frag_' "$NEW_RESULTS" || true)  vs  $(grep -c '_frag_' "$PREV_RESULTS" || true)"
      echo "    first differences:"
      diff "$NEW_RESULTS" "$PREV_RESULTS" | head -20
      FAIL=1
    fi
  fi
fi

hr
if [[ $FAIL -eq 0 ]]; then
  echo "No differences found."
else
  echo "Differences found. If you passed the original --mmf and -s, investigate"
  echo "before publishing: this release should be output-identical on your data."
fi
exit $FAIL
