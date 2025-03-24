# Define the path for the BAM file and the reference name
bam_path = 'bam/fake_cfDNA.bam'
ref_name = 'chr1'

# Define the header for the BAM file
header = {
    'HD': {'VN': '1.0'},
    'SQ': [{'LN': 1000000, 'SN': ref_name}]
}

# Function to generate a random DNA sequence
def random_dna_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

# Create the BAM file with the specified header
with pysam.AlignmentFile(bam_path, 'wb', header=header) as outf:
    for i in range(1000):  # Generate 1000 fragments
        fragment_length = 150 + i % 150  # Vary fragment length between 150 and 300

        # Create the first read of the pair
        a1 = pysam.AlignedSegment()
        a1.query_name = f'read{i}/1'
        a1.query_sequence = random_dna_sequence(100)  # 100 bp read with random sequence
        a1.flag = 99  # Read paired, properly aligned, first in pair
        a1.reference_id = 0
        a1.reference_start = 1000 + i * 200  # Vary start position
        a1.mapping_quality = 60
        a1.cigar = [(0, 100)]  # 100M
        a1.next_reference_id = 0
        a1.next_reference_start = a1.reference_start + fragment_length - 100
        a1.template_length = fragment_length
        a1.query_qualities = [random.randint(20, 40) for _ in range(100)]  # Random quality scores

        # Create the second read of the pair
        a2 = pysam.AlignedSegment()
        a2.query_name = f'read{i}/2'
        a2.query_sequence = random_dna_sequence(100)  # 100 bp read with random sequence
        a2.flag = 147  # Read paired, properly aligned, second in pair
        a2.reference_id = 0
        a2.reference_start = a1.reference_start + fragment_length - 100
        a2.mapping_quality = 60
        a2.cigar = [(0, 100)]  # 100M
        a2.next_reference_id = 0
        a2.next_reference_start = a1.reference_start
        a2.template_length = -fragment_length
        a2.query_qualities = [random.randint(20, 40) for _ in range(100)]  # Random quality scores

        # Write both reads to the BAM file
        outf.write(a1)
        outf.write(a2)

# Index the BAM file
pysam.index(bam_path)
print(f'Generated {bam_path} and its index.')