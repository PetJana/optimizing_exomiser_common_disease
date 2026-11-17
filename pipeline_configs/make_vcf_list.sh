#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Lists the per-sample VCFs of the selected participants.
#
# Usage:
#   ./make_vcf_list.sh <vcf_root_dir> <sample_ids.txt> > samples/<name>_vcfs.txt
#
#   vcf_root_dir    folder searched (recursively) for <ID>.hard-filtered.vcf.gz
#   sample_ids.txt  one VCF sample ID per line (e.g. A113_S1)
#
# Writes one VCF path per line, for run_exomiser.sh. IDs without a VCF are
# reported on stderr.
# ------------------------------------------------------------------------------
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <vcf_root_dir> <sample_ids.txt>" >&2
  exit 1
fi

vcf_root=$1
ids=$2

# All VCFs once, then match on the exact file name (ID followed by the suffix,
# so that e.g. A113_S1 does not also match A113_S12)
all_vcfs=$(find "$vcf_root" -name '*.hard-filtered.vcf.gz')

while IFS= read -r id || [ -n "$id" ]; do
  id=${id%$'\r'}              # tolerate Windows line endings
  [ -z "$id" ] && continue
  match=$(grep -F "/${id}.hard-filtered.vcf.gz" <<< "$all_vcfs" || true)
  if [ -z "$match" ]; then
    echo "No VCF found for $id" >&2
  else
    echo "$match"
  fi
done < "$ids"
