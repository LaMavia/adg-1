#!/usr/bin/env bash

offset=0
len=-1
reads_file=
mapping_file=
dist=
randomise=false
case_name=

function print_usage() {
      cat <<EOM
usege: slice_case.sh [--offset] [--len] [--rand] --input|-i reads_file --output|-o case_name
arguments:
  --input  -i PATH
      Path to the .fasta file with reads.
  --output -o STRING
      Prefix of the file names for the new case: case_name.txt, case_name.fasta.

  --offset    UINT
      Number of skipped cases, non-negative.
  --len       INT
      Max number of cases. Use -1 for all the cases after --offset.
  --rand       
      Randomize reads before slicing.
EOM
}

while [[ $# -gt 0 ]]; do
  case $1 in
    --offset)
      offset=$2
      shift
      shift
      ;;
    --len)
      len=$2
      shift
      shift
      ;;
    --rand)
      randomise=true
      shift
      ;;
    --input | -i)
      reads_file="$(realpath $2)"
      mapping_file="${reads_file%.fasta}.txt"
      dist="$(dirname "$reads_file")"
      shift
      shift
      ;;
    --output | -o)
      case_name="$2"
      shift
      shift
      ;;
    --help)
      print_usage
      exit 0
      ;;
    -* | --*)
      echo "Unknown option $1"
      print_usage
      exit 1
      ;;
    *)
      echo "Unknown positional argument «$1»"
      print_usage
      exit 1
      ;;
  esac
done

if [[ -z "$reads_file" ]]; then
  printf 'Reads file not specified, add --input reads_file_path\n'
  print_usage 
  exit 1
fi

if [[ -z "$case_name" ]]; then
  printf 'Case name not specified, add --output case_name\n'
  print_usage 
  exit 1
fi

if [[ ! -f "$reads_file" ]]; then
  printf 'Failed to find reads file «%s»\n' "$reads_file"
  exit 1
fi

if [[ ! -f "$mapping_file" ]]; then
  printf 'Failed to find mappings file «%s»\n' "$mapping_file"
  exit 1
fi

cd $(dirname $0)

selected_case_pattern=$(if [[ "$randomise" = true ]]; then
    shuf "$mapping_file"
  else
    cat "$mapping_file"
  fi | tail -n "+$(( offset + 1 ))" |\
  if [[ "$len" -lt 0 ]]; then
    cat
  else
    head -n "$len"
  fi | tee "$dist/$case_name.txt" |\
  awk '{ arr[i++] = $1 } END { printf "^> (?:"; for(j in arr) { printf "%s|", arr[j] }; print ")$" }')

cat "$reads_file" | grep -P "$selected_case_pattern" -A 1 --no-group-separator > "$dist/$case_name.fasta"
