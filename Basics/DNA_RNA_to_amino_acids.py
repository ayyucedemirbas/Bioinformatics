def create_codon_table():
    codon_table = {
        # Phenylalanine
        'TTT': 'F', 'TTC': 'F', 'UUU': 'F', 'UUC': 'F',
        # Leucine
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        # Serine
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
        # Tyrosine
        'TAT': 'Y', 'TAC': 'Y', 'UAU': 'Y', 'UAC': 'Y',
        # Cysteine
        'TGT': 'C', 'TGC': 'C', 'UGU': 'C', 'UGC': 'C',
        # Tryptophan
        'TGG': 'W', 'UGG': 'W',
        # Proline
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        # Histidine
        'CAT': 'H', 'CAC': 'H', 'CAU': 'H', 'CAC': 'H',
        # Glutamine
        'CAA': 'Q', 'CAG': 'Q', 'CAA': 'Q', 'CAG': 'Q',
        # Arginine
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        # Isoleucine
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
        # Methionine (Start codon)
        'ATG': 'M', 'AUG': 'M',
        # Threonine
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        # Asparagine
        'AAT': 'N', 'AAC': 'N', 'AAU': 'N', 'AAC': 'N',
        # Lysine
        'AAA': 'K', 'AAG': 'K', 'AAA': 'K', 'AAG': 'K',
        # Valine
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        # Alanine
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        # Aspartic acid
        'GAT': 'D', 'GAC': 'D', 'GAU': 'D', 'GAC': 'D',
        # Glutamic acid
        'GAA': 'E', 'GAG': 'E', 'GAA': 'E', 'GAG': 'E',
        # Glycine
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        # Stop codons
        'TAA': '*', 'TAG': '*', 'TGA': '*', 'UAA': '*', 'UAG': '*', 'UGA': '*'
    }
    return codon_table

def dna_to_rna(dna_sequence):
    return dna_sequence.upper().replace('T', 'U')

def validate_sequence(sequence, seq_type="DNA/RNA"):
    valid_dna = set('ATCG')
    valid_rna = set('AUCG')
    sequence = sequence.upper().replace(' ', '').replace('\n', '')
    
    if seq_type == "DNA":
        valid_chars = valid_dna
    elif seq_type == "RNA":
        valid_chars = valid_rna
    else:
        valid_chars = valid_dna.union(valid_rna)
    
    invalid_chars = set(sequence) - valid_chars
    if invalid_chars:
        raise ValueError(f"Invalid characters found: {invalid_chars}")
    
    return sequence

def translate_sequence(sequence, seq_type="DNA", reading_frame=0):
    sequence = validate_sequence(sequence, seq_type)

    if seq_type.upper() == "DNA":
        rna_sequence = dna_to_rna(sequence)
    else:
        rna_sequence = sequence.upper()

    rna_sequence = rna_sequence[reading_frame:]
    
    codon_table = create_codon_table()
    
    amino_acids = []
    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, 'X')  # X for unknown codon
            amino_acids.append(amino_acid)
    
    return ''.join(amino_acids)

def find_orfs(sequence, seq_type="DNA", min_length=10):
    orfs = []
    
    for frame in range(3):
        translation = translate_sequence(sequence, seq_type, frame)
        
        start_pos = 0
        while True:
            start = translation.find('M', start_pos)
            if start == -1:
                break
            
            stop = translation.find('*', start)
            if stop == -1:
                orf_seq = translation[start:]
                if len(orf_seq) >= min_length:
                    orfs.append({
                        'frame': frame + 1,
                        'start': start * 3 + frame + 1,
                        'end': len(sequence),
                        'length': len(orf_seq),
                        'sequence': orf_seq
                    })
                break
            else:
                orf_seq = translation[start:stop]
                if len(orf_seq) >= min_length:
                    orfs.append({
                        'frame': frame + 1,
                        'start': start * 3 + frame + 1,
                        'end': stop * 3 + frame + 3,
                        'length': len(orf_seq),
                        'sequence': orf_seq
                    })
            
            start_pos = start + 1
    
    return orfs

def main():
    dna_example = "ATGGCATGCAAAGCCTGA"
    rna_example = "AUGUGCAUGCAAAGCCUGA"
    
    print("DNA/RNA to Amino Acid Converter")
    print("=" * 40)
    
    print(f"DNA sequence: {dna_example}")
    dna_translation = translate_sequence(dna_example, "DNA")
    print(f"Translation:  {dna_translation}")
    print(f"Full form:    {' '.join([aa for aa in dna_translation])}")
    
    print(f"RNA sequence: {rna_example}")
    rna_translation = translate_sequence(rna_example, "RNA")
    print(f"Translation:  {rna_translation}")
    
    print(f"All reading frames for DNA sequence:")
    for frame in range(3):
        translation = translate_sequence(dna_example, "DNA", frame)
        print(f"Frame {frame + 1}: {translation}")
    
    # Find ORFs
    longer_sequence = "ATGGCATGCAAAGCCTGATAATGGAGCCCAAATGA"
    print(f"Finding ORFs in: {longer_sequence}")
    orfs = find_orfs(longer_sequence, "DNA", min_length=3)
    
    if orfs:
        print("Open Reading Frames found:")
        for i, orf in enumerate(orfs, 1):
            print(f"ORF {i}: Frame {orf['frame']}, "
                  f"Position {orf['start']}-{orf['end']}, "
                  f"Length {orf['length']} aa")
            print(f"  Sequence: {orf['sequence']}")
    else:
        print("No ORFs found meeting minimum length requirement.")

if __name__ == "__main__":
    main()
