
# Fragment End Analysis Module

This module performs fragment end analysis on cell-free DNA (cfDNA) from BAM files. It extracts key fragmentomic signals relevant to cancer detection.


## Steps
1. Read BAM file and extract fragment lengths, start and end positions.
2. If reference genome is provided, extract k-mers from fragment ends.
3. Calculate motif diversity score (MDS) using normalized Shannon entropy.
4. Normalize motif counts to background and calculate enrichment.
5. Generate summary statistics and save to CSV.
6. Generate and save plots for visualization.

## Features
- Input files: BAM file, reference genome FASTA file (optional)
  - If no reference is provided, end motif distribution will be skipped.
- Main parameters: k-mer length
- Output:
  - `processed.bam`
  - `summary_statistics.csv`
  - `start_motif_distribution.csv` (if reference provided)
  - `end_motif_distribution.csv` (if reference provided)
  - `start_motif_enrichment.csv` (if reference provided)
  - `end_motif_enrichment.csv` (if reference provided)
  - `start_mds_vs_enrichment.png` (if reference provided)
  - `end_mds_vs_enrichment.png` (if reference provided)
  - `motif_diversity_score.txt`
  - `fragment_length_distribution.png`
  - `start_position_distribution.png`
  - `end_position_distribution.png`


## BAM file

For testing purposes, we will use the following publicly available samples from the Roadmap Epigenomic Project.
These were downloaded from nloyfer/wgbs_tools GitHub tutorial folder.
The fastq files were downloaded from GEO, mapped to hg19 using bwa-meth, and sliced to the region chr3:119,527,929-119,531,943 during format conversion (see next section).

SRX | Tissue | Donor
--- | ------ | -----
SRX175350 | Lung cells | STL002
SRX388743 | Pancreas cells | STL002
SRX190161 | Sigmoid colon cells | STL003

### Synthetic data (fake_cfDNA.bam)
The script generate_synthetic.py can be used to generate a synthetic BAM file for testing purposes.

Here is a breakdown of each field:
- Read name: read84/1 - The name of the read, with /1 indicating it is the first read in a pair.
- Flag: 99 - A bitwise flag indicating various properties of the read (e.g., paired, properly aligned, first in pair).
- Reference name: chr1 - The name of the reference sequence (chromosome) to which the read is aligned.
- Position: 17801 - The 1-based leftmost position of the clipped sequence.
- Mapping quality: 60 - The mapping quality score.
- CIGAR string: 100M - The CIGAR string representing the alignment (100 matches).
- Mate reference name: chr1 - The reference name of the mate read.
- Mate position: 17935 - The position of the mate read.
- Template length: 234 - The observed template length.
- Sequence: AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA - The read sequence.
- Quality scores: array('B', [30, 30, 30, ...]) - The quality scores for each base in the read.
- Optional fields: [] - Any additional optional fields (none in this case).


## Reference genome
```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
```

## Requirements
- Python 3.7+
- `pysam`, `pandas`, `seaborn`, `matplotlib`,`numpy`

Install dependencies with:
```bash
pip install -r requirements.txt
```

## Usage

```bash
python fragment_analysis.py --input input.bam --output processed.bam --output_dir results --reference hg38.fa --kmer 3
```


## References

This module implements a minimal fragment end analysis tool inspired by approaches from Snyder et al. (2016). 
For more comprehensive analyses including WPS and TFBS integration, see the Shendure Lab's cfDNA GitHub repository, which this implementation could be extended to emulate.


Calculates a motif diversity score (MDS) using normalized Shannon entropy as described by Jiang et al (2020). This function is generalized for any k instead of just 4-mers.
