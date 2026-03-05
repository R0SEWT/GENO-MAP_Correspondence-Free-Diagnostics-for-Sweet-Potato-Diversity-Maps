#!/usr/bin/env bash
# BLOCKED — ADR-0002: Wild SNP is identical to LowDensity SNP.
#
# The file 01_Report_DSp25-515_SNPs_Filtered_by_reads.csv is an exact copy
# (MD5: ae014cb2a8a17cc215f73b545c571f8b) of the LowDensity SNP file.
# Processing it separately would duplicate 62,736 markers × 635 samples.
#
# Use run_lowdensity_snp.sh instead.
# See: docs/addr/0002-duplicado-lowdensity-wild-snp.md

set -euo pipefail

echo "ERROR: Wild SNP es idéntico a LowDensity SNP (ADR-0002)."
echo "       Usa run_lowdensity_snp.sh en su lugar."
echo "       Detalle: docs/addr/0002-duplicado-lowdensity-wild-snp.md"
exit 1
