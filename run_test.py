import os
import glob
import subprocess

bam_dir = 'cfdx/bam'
results_dir = 'cfdx/results'

os.makedirs(results_dir, exist_ok=True)

for bam_file in glob.glob(os.path.join(bam_dir, '*.bam')):
    base_name = os.path.basename(bam_file).replace('.bam', '')
    bam_results_dir = os.path.join(results_dir, base_name)
    os.makedirs(bam_results_dir, exist_ok=True)

    output_bam = os.path.join(bam_results_dir, f'{base_name}.processed.bam')
    reference_fasta = 'cfdx/hg19.fa'

    subprocess.run(['python', 'cfdx/fragment_analysis.py', '--input', bam_file, '--output', output_bam, '--output_dir', bam_results_dir, '--reference', reference_fasta, '--kmer', '3'])