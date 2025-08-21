import random
from itertools import product

def create_codon_table():
    codon_table = {
        'F': ['TTT', 'TTC'],  # Phenylalanine
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # Leucine
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Serine
        'Y': ['TAT', 'TAC'],  # Tyrosine
        'C': ['TGT', 'TGC'],  # Cysteine
        'W': ['TGG'],  # Tryptophan
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],  # Proline
        'H': ['CAT', 'CAC'],  # Histidine
        'Q': ['CAA', 'CAG'],  # Glutamine
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arginine
        'I': ['ATT', 'ATC', 'ATA'],  # Isoleucine
        'M': ['ATG'],  # Methionine (Start)
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],  # Threonine
        'N': ['AAT', 'AAC'],  # Asparagine
        'K': ['AAA', 'AAG'],  # Lysine
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],  # Valine
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],  # Alanine
        'D': ['GAT', 'GAC'],  # Aspartic acid
        'E': ['GAA', 'GAG'],  # Glutamic acid
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],  # Glycine
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    }
    return codon_table

def validate_amino_acid_sequence(sequence):
    valid_aa = set('FLSYCWPHQRIMTNKVADEG*')
    sequence = sequence.upper().replace(' ', '').replace('\n', '')
    
    invalid_chars = set(sequence) - valid_aa
    if invalid_chars:
        raise ValueError(f"Invalid amino acids found: {invalid_chars}")
    
    return sequence

def amino_acid_to_dna(aa_sequence, method="random"):
    aa_sequence = validate_amino_acid_sequence(aa_sequence)
    codon_table = create_codon_table()
    
    dna_codons = []
    for aa in aa_sequence:
        if aa not in codon_table:
            raise ValueError(f"Unknown amino acid: {aa}")
        
        possible_codons = codon_table[aa]
        
        if method == "random":
            chosen_codon = random.choice(possible_codons)
        elif method == "first":
            chosen_codon = possible_codons[0]
        elif method == "optimal":
            # Use most frequent codons in E. coli
            optimal_codons = {
                'F': 'TTT', 'L': 'CTG', 'S': 'TCG', 'Y': 'TAT', 'C': 'TGC',
                'W': 'TGG', 'P': 'CCG', 'H': 'CAT', 'Q': 'CAG', 'R': 'CGT',
                'I': 'ATT', 'M': 'ATG', 'T': 'ACC', 'N': 'AAC', 'K': 'AAA',
                'V': 'GTG', 'A': 'GCG', 'D': 'GAT', 'E': 'GAA', 'G': 'GGT',
                '*': 'TAA'
            }
            chosen_codon = optimal_codons.get(aa, possible_codons[0])
        
        dna_codons.append(chosen_codon)
    
    return ''.join(dna_codons)

def amino_acid_to_rna(aa_sequence, method="random"):
    dna_seq = amino_acid_to_dna(aa_sequence, method)
    return dna_seq.replace('T', 'U')

def get_all_possible_sequences(aa_sequence, seq_type="DNA", max_combinations=1000):
    aa_sequence = validate_amino_acid_sequence(aa_sequence)
    codon_table = create_codon_table()
    
    codon_options = []
    for aa in aa_sequence:
        codon_options.append(codon_table[aa])
    
    total_combinations = 1
    for options in codon_options:
        total_combinations *= len(options)
    
    print(f"Total possible combinations: {total_combinations}")
    
    if total_combinations > max_combinations:
        print(f"Limiting to first {max_combinations} combinations")
    
    sequences = []
    for i, combination in enumerate(product(*codon_options)):
        if i >= max_combinations:
            break
        
        seq = ''.join(combination)
        if seq_type.upper() == "RNA":
            seq = seq.replace('T', 'U')
        
        sequences.append(seq)
    
    return sequences

def calculate_degeneracy(aa_sequence):
    aa_sequence = validate_amino_acid_sequence(aa_sequence)
    codon_table = create_codon_table()
    
    degeneracies = []
    for i, aa in enumerate(aa_sequence):
        num_codons = len(codon_table[aa])
        degeneracies.append({
            'position': i + 1,
            'amino_acid': aa,
            'possible_codons': num_codons,
            'codons': codon_table[aa]
        })
    
    return degeneracies

def translate_with_constraints(aa_sequence, constraints=None):
    aa_sequence = validate_amino_acid_sequence(aa_sequence)
    codon_table = create_codon_table()
    constraints = constraints or {}
    
    dna_codons = []
    for i, aa in enumerate(aa_sequence):
        if i in constraints:
            preferred_codon = constraints[i].upper()
            if preferred_codon in codon_table[aa]:
                dna_codons.append(preferred_codon)
            else:
                print(f"Warning: {preferred_codon} not valid for {aa} at position {i}")
                dna_codons.append(random.choice(codon_table[aa]))
        else:
            dna_codons.append(random.choice(codon_table[aa]))
    
    return ''.join(dna_codons)

def main():
    print("Amino Acid to DNA/RNA Converter")
    print("=" * 40)
    
    aa_seq = "MACK*"  # Met-Ala-Cys-Lys-Stop
    
    print(f"Amino acid sequence: {aa_seq}")
    
    print("\nDNA translations:")
    for method in ["random", "first", "optimal"]:
        dna_seq = amino_acid_to_dna(aa_seq, method)
        print(f"{method.capitalize():8}: {dna_seq}")
    
    rna_seq = amino_acid_to_rna(aa_seq, "random")
    print(f"\nRNA (random): {rna_seq}")
    
    print(f"\nCodon degeneracy analysis:")
    degeneracies = calculate_degeneracy(aa_seq)
    for deg in degeneracies:
        print(f"Position {deg['position']} ({deg['amino_acid']}): "
              f"{deg['possible_codons']} possible codons")
        print(f"  Options: {', '.join(deg['codons'])}")
    
    print(f"\nFirst 10 possible DNA sequences:")
    possible_seqs = get_all_possible_sequences(aa_seq, "DNA", 10)
    for i, seq in enumerate(possible_seqs, 1):
        print(f"{i:2}: {seq}")
    
    constraints = {0: 'ATG', 4: 'TGA'}
    constrained_seq = translate_with_constraints(aa_seq, constraints)
    print(f"\nWith constraints (pos 0=ATG, pos 4=TGA): {constrained_seq}")

if __name__ == "__main__":
    main()
