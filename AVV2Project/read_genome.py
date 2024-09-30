# Structure: The AAV2 genome contains two main genes:
# Rep (Replication): Involved in the replication of the viral genome.
# Cap (Capsid): Encodes the proteins that form the viral capsid (the protein shell of the virus).

# Define the path to your FASTA file
fasta_file = 'avv2.fna'


# Function to read a FASTA file and return the sequence
def read_fasta(file):
    with open(file, 'r') as f:
        # Skip the header line and join the sequence lines
        seq = ''.join(line.strip() for line in f if not line.startswith('>'))
    return seq


# Load the genome sequence
avv2_genome = read_fasta(fasta_file)

# Define lengths for ITRs
ITR_length = 145

# Extract ITRs
itr_1 = avv2_genome[:ITR_length]  # First ITR
itr_2 = avv2_genome[-ITR_length:]  # Second ITR


# Find all genes that start with ATG and end with TAA, TAG, or TGA
def find_genes_with_start_stop(genome):
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    genes = []  # List to store found genes

    # Start searching from the beginning of the genome
    start_pos = 0
    while True:
        # Find the next start codon
        start_pos = genome.find(start_codon, start_pos)
        if start_pos == -1:
            break  # No more start codons found

        # Now find the next stop codon after the start codon
        stop_pos = -1
        for codon in stop_codons:
            pos = genome.find(codon, start_pos + 3)  # Look for stop codon after the start codon
            if pos != -1:
                if stop_pos == -1 or pos < stop_pos:  # Find the earliest stop codon
                    stop_pos = pos

        if stop_pos != -1:
            # Extract the gene sequence from start to stop codon
            gene_sequence = genome[start_pos:stop_pos + 3]  # Include the stop codon
            genes.append(gene_sequence)  # Append to the list of genes

        # Move past the current start codon for the next search
        start_pos += 3

    return genes


# Find genes in the AAV2 genome
genes = find_genes_with_start_stop(avv2_genome)

# Print the results
print("First ITR: ", itr_1)
print("Second ITR: ", itr_2)

# Print found genes
for i, gene in enumerate(genes, start=1):
    print(f"Gene {i}: ", gene)

# Optionally, write the results to a file
with open('aav2_genome_results.txt', 'w') as output_file:
    output_file.write("First ITR: {}\n".format(itr_1))
    output_file.write("Second ITR: {}\n".format(itr_2))
    for i, gene in enumerate(genes, start=1):
        output_file.write(f"Gene {i}: {gene}\n")
