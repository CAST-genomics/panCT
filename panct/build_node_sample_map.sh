#!/usr/bin/env sh

# Extract walks (W lines) from a GFA file, flip them so that they can be queried by node instead of sample, and write the results to a TSV file

# arg1: GFA file from which to extract walks
# arg2 (optional): Path to .walk tab separated file to which to write output. Defaults to stdout

# ex: ./build_node_sample_map.sh tests/data/basic.gfa basic.walk


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
}' | sort -t '	' -k1,1n -k2,2 | awk -F '	' -v OFS='	' '
current==$1 { line = line OFS $2; next; }
{ print line; current=$1; line=$0; }
END { print line }' | awk -F '	' -v OFS='	' 'NF {
  print "" FS $0
}' | {
  if [ "${2##*.}" = "gz" ]; then
    bgzip
  else
    cat
  fi
} > "${2:-/dev/stdout}"
