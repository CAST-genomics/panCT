#!/usr/bin/env sh

{ if [ "${1##*.}" = "gz" ]; then
  gzip -dc "$1"
else
  cat "$1"
fi } | grep -E '^W' | cut -f 2,7 | sed 's/	//;s/</	/g;s/>/	/g' | awk -F '	' -v OFS='	' '{
  for (i = 1; i <= NF; i++) {
    if (i == 1) {
      key = $i;
    } else {
      print $i, key;
    }
  }
}' | sort -u -t '	' -k1,1n -k2,2 | awk -F '	' -v OFS='	' '
current==$1 { line = line OFS $2; next; }
{ print line; current=$1; line=$0; }
END { print line }' | awk -F '	' -v OFS='	' 'NF {
  print "" FS $1 FS $0
}' > "$2"
