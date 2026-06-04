#!/bin/bash
set -o errexit
set -o pipefail

find src/vstolib test -name '*.py' \
  | xargs pylint \
  --errors-only \
  --disable=unsubscriptable-object,not-an-iterable,no-member

echo 'Passes pylint check'