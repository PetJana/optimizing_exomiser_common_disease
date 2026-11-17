#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Runs one Exomiser analysis for a list of samples, using Exomiser's batch mode.
#
# Usage (from the Exomiser install directory, so that relative outputDirectory
# paths in the configs resolve to <install_dir>/results/<run>):
#   <repo>/pipeline_configs/run_exomiser.sh <analysis.yml> <vcf_list.txt>
#
#   analysis.yml   one of pipeline_configs/analysis/*.yml
#   vcf_list.txt   one VCF path per line (see make_vcf_list.sh)
#
# Environment variables (optional):
#   EXOMISER_JAR   Exomiser CLI jar   (default: exomiser-cli-14.0.0.jar)
#   JAVA_MEM       Java heap size     (default: 4g)
# ------------------------------------------------------------------------------
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <analysis.yml> <vcf_list.txt>" >&2
  exit 1
fi

analysis=$(realpath "$1")
vcf_list=$(realpath "$2")
jar=${EXOMISER_JAR:-exomiser-cli-14.0.0.jar}
mem=${JAVA_MEM:-4g}

# One line of Exomiser arguments per sample
batch_file=$(mktemp)
trap 'rm -f "$batch_file"' EXIT

while IFS= read -r vcf || [ -n "$vcf" ]; do
  vcf=${vcf%$'\r'}
  [ -z "$vcf" ] && continue
  echo "--analysis $analysis --vcf $vcf --assembly hg19"
done < "$vcf_list" > "$batch_file"

echo "Running $(basename "$analysis") on $(wc -l < "$batch_file") samples"
java -Xmx"$mem" -jar "$jar" --batch "$batch_file"
