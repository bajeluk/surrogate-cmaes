#!/bin/sh

# test [bbobexp_f[FUN].info] files whether the 3-rd line ends with 40:XXX|YYYYY
#
# expects directory(-ies) as the input parameter(s)

FILE='bbobexp_f*.info'
IS40_REGEXP='3s/.*\( 40:[0-9]\{2,5\}|[^ ]\+\)$/\1/p'

for f in $*; do
  is40last=`sed -n "$IS40_REGEXP" $f/$FILE`
  if [ -z "$is40last" ]; then
    echo $f ' <-- !!!'
  fi
done
