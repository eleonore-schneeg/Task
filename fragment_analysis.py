import os
import pysam
import pandas as pd
from collections import Counter
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def filter_read(read, min_quality=20, min_length=30):
    """Filter reads that do not meet quality criteria."""
    if read.query_qualities is None:
        return False  # No quality scores to filter

    if len(read.query_sequence) < min_length:
        return False
    if min(read.query_qualities) < min_quality:
        return False
    return True

def motif_diversity_score(motif_counts: dict) -> float:
    total = sum(motif_counts.values())
    if total == 0:
        return 0.0
    probs = np.array([count / total for count in motif_counts.values()])
    entropy = -np.sum(probs * np.log2(probs))
    max_entropy = np.log2(len(motif_counts)) if motif_counts else 1
    return entropy / max_entropy if max_entropy > 0 else 0.0

def normalize_to_background(motif_counts: dict, kmer_size: int) -> pd.DataFrame:
    all_kmers = ["".join(p) for p in __import__('itertools').product("ACGT", repeat=kmer_size)]
    observed = pd.Series(motif_counts, index=all_kmers).fillna(0)
    expected = pd.Series(1.0 / len(all_kmers), index=all_kmers) * observed.sum()
    enrichment = observed / expected
    return pd.DataFrame({"Motif": enrichment.index, "Observed": observed.values, "Expected": expected.values, "Enrichment": enrichment.values})

def plot_mds_vs_enrichment(enrichment_df: pd.DataFrame, title: str, output_file: str):
    enrichment_df = enrichment_df.copy()
    enrichment_df["Log2Enrichment"] = np.log2(enrichment_df["Enrichment"].replace(0, np.nan))
    enrichment_df.dropna(subset=["Log2Enrichment"], inplace=True)
    enrichment_df = enrichment_df.sort_values("Log2Enrichment", ascending=False).head(5)
    enrichment_df["Rank"] = range(1, len(enrichment_df) + 1)

    plt.figure(figsize=(12, 8))
    ax = sns.scatterplot(x="Rank", y="Log2Enrichment", data=enrichment_df)

    for i, row in enrichment_df.iterrows():
        ax.text(row["Rank"], row["Log2Enrichment"], row["Motif"], fontsize=9, ha='center', va='bottom')

    plt.title(title)
    plt.xlabel("Motif Rank")
    plt.ylabel("log2(Enrichment)")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def analyze_bam(input_bam, output_bam, output_dir, reference_fasta=None, kmer_size=3, min_quality=20, min_length=30):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    frag_lengths = []
    start_pos = []
    end_pos = []
    start_motifs = Counter()
    end_motifs = Counter()

    if reference_fasta:
        ref = pysam.FastaFile(reference_fasta)

    for read in bam_in.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        if not filter_read(read, min_quality=min_quality, min_length=min_length):
            continue

        start = read.reference_start
        end = read.reference_end
        fragment_length = end - start

        frag_lengths.append(fragment_length)
        start_pos.append(start)
        end_pos.append(end)

        if reference_fasta:
            try:
                if read.is_reverse:
                    start_motif = ref.fetch(read.reference_name, end - kmer_size, end)
                    start_motif = str(Seq(start_motif).reverse_complement())
                    end_motif = ref.fetch(read.reference_name, start, start + kmer_size)
                    end_motif = str(Seq(end_motif).reverse_complement())
                else:
                    start_motif = ref.fetch(read.reference_name, start, start + kmer_size)
                    end_motif = ref.fetch(read.reference_name, end - kmer_size, end)

                start_motifs[start_motif.upper()] += 1
                end_motifs[end_motif.upper()] += 1
            except:
                continue

        bam_out.write(read)

    bam_in.close()
    bam_out.close()

    os.makedirs(output_dir, exist_ok=True)

    summary_df = pd.DataFrame({"FragmentLength": frag_lengths, "Start": start_pos, "End": end_pos})
    summary_df.to_csv(os.path.join(output_dir, "summary_statistics.csv"), index=False)

    if start_motifs:
        pd.DataFrame(start_motifs.items(), columns=["Motif", "Count"]).to_csv(os.path.join(output_dir, "start_motif_distribution.csv"), index=False)
        enriched_start = normalize_to_background(start_motifs, kmer_size)
        enriched_start.to_csv(os.path.join(output_dir, "start_motif_enrichment.csv"), index=False)
        plot_mds_vs_enrichment(enriched_start, "5' End: log2 Enrichment by Motif Rank", os.path.join(output_dir, "start_mds_vs_enrichment.png"))

    if end_motifs:
        pd.DataFrame(end_motifs.items(), columns=["Motif", "Count"]).to_csv(os.path.join(output_dir, "end_motif_distribution.csv"), index=False)
        enriched_end = normalize_to_background(end_motifs, kmer_size)
        enriched_end.to_csv(os.path.join(output_dir, "end_motif_enrichment.csv"), index=False)
        plot_mds_vs_enrichment(enriched_end, "3' End: log2 Enrichment by Motif Rank", os.path.join(output_dir, "end_mds_vs_enrichment.png"))

    mds_start = motif_diversity_score(start_motifs)
    mds_end = motif_diversity_score(end_motifs)
    with open(os.path.join(output_dir, "motif_diversity_score.txt"), "w") as f:
        f.write(f"5' End Motif Diversity Score: {mds_start:.4f}\n")
        f.write(f"3' End Motif Diversity Score: {mds_end:.4f}\n")

    plt.figure()
    sns.histplot(frag_lengths, bins=100, kde=True)
    plt.title("Fragment Length Distribution")
    plt.xlabel("Fragment Length")
    plt.ylabel("Count")
    plt.savefig(os.path.join(output_dir, "fragment_length_distribution.png"))

    plt.figure()
    sns.histplot(start_pos, bins=100, kde=False)
    plt.title("Start Position Distribution")
    plt.xlabel("Genomic Start Position")
    plt.ylabel("Count")
    plt.savefig(os.path.join(output_dir, "start_position_distribution.png"))

    plt.figure()
    sns.histplot(end_pos, bins=100, kde=False)
    plt.title("End Position Distribution")
    plt.xlabel("Genomic End Position")
    plt.ylabel("Count")
    plt.savefig(os.path.join(output_dir, "end_position_distribution.png"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fragment End Analysis on cfDNA BAM files.")
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output processed BAM file")
    parser.add_argument("--output_dir", required=True, help="Output directory for results")
    parser.add_argument("--reference", required=False, help="Reference FASTA for motif extraction")
    parser.add_argument("--kmer", type=int, default=3, help="Length of k-mer to extract from fragment ends")
    parser.add_argument("--min_quality", type=int, default=20, help="Minimum quality score for filtering reads")
    parser.add_argument("--min_length", type=int, default=20, help="Minimum length for filtering reads")

    args = parser.parse_args()
    analyze_bam(args.input, args.output, args.output_dir, args.reference, args.kmer, args.min_quality)