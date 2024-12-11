#!/usr/bin/env sh

zcat "$1" | grep -E '^W' | cut -f 2,7 | sed 's/\t//;s/</\t/g;s/>/\t/g' | awk -F $'\t' -v 'OFS=\t' '{
  for (i = 1; i <= NF; i++) {
    if (i == 1) {
      key = $i;
    } else {
      print($i, key);
    }
  }
}' | sort -u -t $'\t' -k1,1n -k2,2 | awk -F $'\t' -v 'OFS=\t' '
current==$1 { line = line OFS $2; next; }
{ print line; current=$1; line=$0; }
END { print line }' | awk -F $'\t' -v 'OFS=\t' 'NF {
  print "" FS $1 FS $0
}' > "$2"
