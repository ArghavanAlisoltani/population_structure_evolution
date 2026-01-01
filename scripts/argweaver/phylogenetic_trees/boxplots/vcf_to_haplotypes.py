import sys

def extract_haplotypes(vcf_path, target_scaffold, target_positions):
    """
    Parses a VCF file and extracts alleles for specific positions 
    to create two haplotypes per sample.
    """
    samples = []
    # Dictionary to store genotypes: {sample_index: [hap1_list, hap2_list]}
    haplotypes = {}
    
    try:
        with open(vcf_path, 'r') as f:
            for line in f:
                # Skip irrelevant header lines
                if line.startswith("##"):
                    continue
                
                # Parse the column header line to get sample names
                if line.startswith("#CHROM"):
                    header = line.strip().split('\t')
                    samples = header[9:] # Samples start at the 10th column
                    for i in range(len(samples)):
                        haplotypes[i] = [[], []]
                    continue
                
                # Parse data lines
                cols = line.strip().split('\t')
                scaffold = cols[0]
                pos = int(cols[1])
                ref = cols[3]
                alt_alleles = cols[4].split(',')
                all_alleles = [ref] + alt_alleles
                
                # Check if this row matches our criteria
                if scaffold == target_scaffold and pos in target_positions:
                    format_col = cols[8].split(':')
                    gt_index = format_col.index('GT') if 'GT' in format_col else 0
                    
                    sample_data = cols[9:]
                    for i, data in enumerate(sample_data):
                        # Extract GT field (e.g., "0|1" or "0/1")
                        gt_field = data.split(':')[gt_index]
                        
                        # Handle phased (|) or unphased (/) separators
                        sep = '|' if '|' in gt_field else '/'
                        alleles = gt_field.split(sep)
                        
                        # Map index to actual base (handle missing data '.' as 'N')
                        base1 = all_alleles[int(alleles[0])] if alleles[0] != '.' else 'N'
                        base2 = all_alleles[int(alleles[1])] if alleles[1] != '.' else 'N'
                        
                        haplotypes[i][0].append(base1)
                        haplotypes[i][1].append(base2)

        # Output in FASTA format
        for i, sample_name in enumerate(samples):
            # Haplotype 1
            print(f">{sample_name}_1")
            print("".join(haplotypes[i][0]))
            # Haplotype 2
            print(f">{sample_name}_2")
            print("".join(haplotypes[i][1]))

    except FileNotFoundError:
        print(f"Error: File {vcf_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # --- CONFIGURATION ---
    # Update these values to match your requirements
    input_vcf = "vcf.header.txt"
    scaffold = "scaffold_17"
    positions = [21197, 21238, 21240, 21245, 21297, 21300, 21302] 
    # (The example positions are from your provided file content)
    
    extract_haplotypes(input_vcf, scaffold, positions)
