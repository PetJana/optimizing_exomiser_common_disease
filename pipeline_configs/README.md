# Exomiser configuration

Exomiser analysis settings and helper scripts used to produce the per-participant Exomiser results that the R pipeline in the parent folder reads.

## Contents

```
pipeline_configs/
├── analysis/                              # one Exomiser analysis file per run
│   ├── SKSC_3_MAF1per_herIncluded_onlyKS.yml
│   ├── SKSC_3_MAF1per_herIncluded.yml
│   ├── SKSC_3_MAF1per.yml
│   └── SKSC_3_MAF1per_nonHIT.yml
├── hpo_terms_SKSC.tsv                     # the kidney stone-related HPO terms, with labels
├── make_vcf_list.sh                       # sample IDs -> list of VCF paths
├── run_exomiser.sh                        # runs one analysis on a list of VCFs (batch mode)
└── samples/                               # (not tracked) sample ID and VCF lists
```

## Runs

Each analysis file writes to a folder with the same name as the run. These are the folders the R pipeline expects (`exomiser_dir` in `00_setup.R`).

| Run | HPO terms | Inheritance-mode weighting | Samples | Used for |
|---|---|---|---|---|
| `SKSC_3_MAF1per_herIncluded_onlyKS` | nephrolithiasis (HP:0000787) only | yes (`omimPrioritiser`) | benchmark participants | benchmarking — "KSD, default settings" |
| `SKSC_3_MAF1per_herIncluded` | kidney stone-related terms | yes (`omimPrioritiser`) | benchmark participants | benchmarking — "KSD+, default settings" |
| `SKSC_3_MAF1per` | kidney stone-related terms | no | benchmark participants | benchmarking — "KSD+, any inh."; ESTABLISHED genes |
| `SKSC_3_MAF1per_nonHIT` | kidney stone-related terms | no | all participants | NON-ESTABLISHED genes (candidate gene nomination) |

Settings shared by all runs:

- genome assembly GRCh37 (hg19), PASS variants only
- maximum minor allele frequency 1% for dominant and homozygous recessive models (2% for compound heterozygous, 0.2% mitochondrial), using UK10K and gnomAD exomes/genomes (AFR, AMR, EAS, NFE, SAS)
- pathogenicity sources REVEL, MVP, SpliceAI and AlphaMissense
- UTR, intronic, non-coding, upstream/downstream, intergenic and regulatory variants removed
- gene prioritisation with hiPHIVE

The kidney stone-related HPO terms are listed in `hpo_terms_SKSC.tsv`. 

## Running Exomiser

Requires Exomiser v14.0.0 (CLI) with the 2406 data release (hg19 and phenotype data), configured in Exomiser's `application.properties` (including the SpliceAI and AlphaMissense sources). See the [Exomiser documentation](https://exomiser.readthedocs.io/).

1. Make a list of sample IDs for each run (one VCF sample ID per line), e.g. `samples/benchmark_ids.txt` and `samples/all_ids.txt`.
2. Turn each into a list of VCF paths:

   ```bash
   ./make_vcf_list.sh /path/to/vcfs samples/benchmark_ids.txt > samples/benchmark_vcfs.txt
   ./make_vcf_list.sh /path/to/vcfs samples/all_ids.txt       > samples/all_vcfs.txt
   ```

3. From the Exomiser install directory, run each analysis:

   ```bash
   cd /path/to/exomiser-cli-14.0.0
   CFG=/path/to/repo/pipeline_configs
   $CFG/run_exomiser.sh $CFG/analysis/SKSC_3_MAF1per_herIncluded_onlyKS.yml $CFG/samples/benchmark_vcfs.txt
   $CFG/run_exomiser.sh $CFG/analysis/SKSC_3_MAF1per_herIncluded.yml        $CFG/samples/benchmark_vcfs.txt
   $CFG/run_exomiser.sh $CFG/analysis/SKSC_3_MAF1per.yml                    $CFG/samples/benchmark_vcfs.txt
   $CFG/run_exomiser.sh $CFG/analysis/SKSC_3_MAF1per_nonHIT.yml             $CFG/samples/all_vcfs.txt
   ```

Results are written to `results/<run>/` in the Exomiser install directory. Point `exomiser_dir` in `00_setup.R` to the folder that contains these run folders.

The sample and VCF lists contain participant IDs and are not part of the repository.
