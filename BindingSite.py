#!/usr/bin/env python3

import os
import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

class DNAShapePredictor:
    """Predicts DNA shape features"""
    def __init__(self):
        self.shape_params = {
            'AA': {'mgw': 4.0, 'prop_tw': -14.0, 'roll': 0.6},
            'AT': {'mgw': 4.2, 'prop_tw': -13.3, 'roll': 1.1},
            'AC': {'mgw': 4.1, 'prop_tw': -14.5, 'roll': 0.9},
            'AG': {'mgw': 4.0, 'prop_tw': -14.0, 'roll': 1.0},
            'TA': {'mgw': 4.7, 'prop_tw': -11.8, 'roll': 3.3},
            'TT': {'mgw': 4.0, 'prop_tw': -14.0, 'roll': 0.6},
            'TC': {'mgw': 4.1, 'prop_tw': -14.8, 'roll': 0.7},
            'TG': {'mgw': 4.2, 'prop_tw': -14.3, 'roll': 1.2},
            'CA': {'mgw': 4.1, 'prop_tw': -14.8, 'roll': 0.7},
            'CT': {'mgw': 4.1, 'prop_tw': -14.5, 'roll': 0.9},
            'CC': {'mgw': 4.2, 'prop_tw': -15.0, 'roll': 0.8},
            'CG': {'mgw': 4.3, 'prop_tw': -14.8, 'roll': 1.3},
            'GA': {'mgw': 4.0, 'prop_tw': -14.0, 'roll': 1.0},
            'GT': {'mgw': 4.2, 'prop_tw': -14.3, 'roll': 1.2},
            'GC': {'mgw': 4.3, 'prop_tw': -14.8, 'roll': 1.3},
            'GG': {'mgw': 4.2, 'prop_tw': -14.5, 'roll': 1.1}
        }

    def predict_shape(self, sequence):
        try:
            features = {'minor_groove_width': [], 'propeller_twist': [], 'roll': []}
            if len(sequence) < 2:
                return {k: [0] for k in features.keys()}

            for i in range(len(sequence)-1):
                dinuc = sequence[i:i+2].upper()
                if dinuc in self.shape_params:
                    params = self.shape_params[dinuc]
                    for key in features:
                        features[key].append(params[key.split('_')[0] if '_' in key else key])

            if not any(features.values()):
                return {k: [0] for k in features.keys()}
            return features
        except Exception as e:
            print(f"Error in predict_shape: {e}")
            return {k: [0] for k in features.keys()}

class SinRBindingPredictor:
    """Predicts SinR binding sites with p-values"""
    def __init__(self):
        self.known_motifs = [
            "GTTCTCT",
            "AGAAGAC",
            "GTTNNNNNNNNAAC",
            "CACGAAAT",
            "TGAAAT"
        ]
        self.shape_predictor = DNAShapePredictor()
        self.null_distribution = None
        self.null_params = None
        self.null_distribution_generated = False

    def generate_null_distribution(self, background_sequences, num_samples=10000):
        """Generate null distribution of binding scores"""
        null_scores = []
        for _ in tqdm(range(num_samples), desc="Generating null distribution"):
            seq = np.random.choice(background_sequences)
            score = self.predict_binding_affinity(seq)
            null_scores.append(score)
        
        self.null_distribution = np.array(null_scores)
        self.null_params = stats.genextreme.fit(self.null_distribution)
        self.null_distribution_generated = True
        return self.null_distribution

    def calculate_pvalue(self, score):
        """Calculate p-value for a given binding score"""
        if not self.null_distribution_generated:
            return np.nan
        return stats.genextreme.sf(score, *self.null_params)

    def _calculate_palindrome_score(self, sequence):
        try:
            if not sequence:
                return 0
            rev_comp = str(Seq(sequence).reverse_complement())
            matches = sum(a == b for a, b in zip(sequence, rev_comp))
            return matches / len(sequence)
        except Exception:
            return 0

    def _calculate_conservation_score(self, sequence):
        try:
            if not sequence:
                return 0
            scores = []
            for motif in self.known_motifs:
                if len(sequence) < len(motif):
                    continue
                if 'N' in motif:
                    non_n_positions = [i for i, c in enumerate(motif) if c != 'N']
                    motif_bases = [motif[i] for i in non_n_positions]
                    seq_bases = [sequence[i] for i in non_n_positions if i < len(sequence)]
                    if seq_bases:
                        score = sum(a == b for a, b in zip(motif_bases, seq_bases)) / len(seq_bases)
                        scores.append(score)
                else:
                    matches = sum(a == b for a, b in zip(sequence, motif))
                    scores.append(matches / len(motif))
            return max(scores) if scores else 0
        except Exception:
            return 0

    def calculate_shape_score(self, sequence):
        try:
            shape_features = self.shape_predictor.predict_shape(sequence)
            if not any(shape_features.values()):
                return 0

            avg_mgw = np.mean(shape_features['minor_groove_width'])
            avg_ptw = np.mean(shape_features['propeller_twist'])
            avg_roll = np.mean(shape_features['roll'])

            shape_score = (
                0.4 * (avg_mgw / 5.0) +
                0.3 * (abs(avg_ptw) / 15.0) +
                0.3 * (avg_roll / 3.5)
            )
            return shape_score
        except Exception:
            return 0

    def predict_binding_affinity(self, sequence):
        try:
            if not sequence or len(sequence) < 6:
                return 0

            at_content = (sequence.count('A') + sequence.count('T')) / len(sequence)
            palindrome_score = self._calculate_palindrome_score(sequence)
            conservation_score = self._calculate_conservation_score(sequence)
            shape_score = self.calculate_shape_score(sequence)

            score = (
                0.25 * at_content +
                0.25 * palindrome_score +
                0.25 * conservation_score +
                0.25 * shape_score
            )
            return score
        except Exception:
            return 0

    def predict_binding_affinity_with_pvalue(self, sequence):
        score = self.predict_binding_affinity(sequence)
        pvalue = self.calculate_pvalue(score)
        return score, pvalue

def get_user_input(prompt, is_file=False):
    while True:
        value = input(prompt).strip()
        if not value:
            print("Input cannot be empty. Please try again.")
            continue
        if is_file:
            if not os.path.exists(value):
                print(f"File '{value}' does not exist. Please try again.")
                continue
            try:
                with open(value) as f:
                    first_line = f.readline()
                    if not first_line.startswith('>'):
                        print("File doesn't appear to be in FASTA format. Please provide a valid FASTA file.")
                        continue
            except Exception as e:
                print(f"Error reading file: {e}")
                continue
        return value

def read_fasta_file(file_path):
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            return record.id, str(record.seq)
        raise ValueError("No sequences found in FASTA file")
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)

def process_sequence_chunk(chunk, predictor):
    try:
        sequence = chunk['sequence']
        window_size = 20
        step_size = 10
        results = []
        seen = set()

        for i in range(0, len(sequence) - window_size + 1, step_size):
            window = sequence[i:i + window_size]
            position = chunk['start'] + i
            key = (position, window)

            if key in seen:
                continue
            seen.add(key)

            score, pvalue = predictor.predict_binding_affinity_with_pvalue(window)
            results.append({
                'position': position,
                'sequence': window,
                'score': score,
                'pvalue': pvalue,
                'contig': chunk['contig']
            })

        return results
    except Exception as e:
        print(f"Error processing chunk: {e}")
        return []

def main():
    RANDOM_SEED = 42
    np.random.seed(RANDOM_SEED)

    print("\n=== SinR Binding Site Analysis Tool ===\n")

    print("Step 1: Protein Information")
    print("-" * 30)
    protein_file = get_user_input("Enter path to protein FASTA file: ", is_file=True)
    protein_name, protein_sequence = read_fasta_file(protein_file)
    print(f"\nLoaded protein: {protein_name}")
    print(f"Sequence length: {len(protein_sequence)} amino acids")

    print("\nStep 2: Genome Information")
    print("-" * 30)
    genome_file = get_user_input("Enter path to genome FASTA file: ", is_file=True)
    genome_name = get_user_input("Enter genome name: ")

    output_dir = f"{protein_name}_{genome_name}_analysis"
    os.makedirs(output_dir, exist_ok=True)

    with open(os.path.join(output_dir, "seed_used.txt"), "w") as f:
        f.write(f"Random seed used: {RANDOM_SEED}\n")

    print("\nStep 3: Analysis Parameters")
    print("-" * 30)
    chunk_size = 1000
    step_size = 500
    num_processes = os.cpu_count() or 4
    print(f"Using {num_processes} CPU cores for parallel processing")

    print("\nStep 4: Processing Genome")
    print("-" * 30)
    chunks = []
    genome_size = 0

    print("Reading genome file...")
    for record in SeqIO.parse(genome_file, "fasta"):
        seq = str(record.seq)
        genome_size += len(seq)
        for i in range(0, len(seq) - chunk_size + 1, step_size):
            chunks.append({
                'contig': record.id,
                'start': i,
                'sequence': seq[i:i + chunk_size]
            })

    print(f"\nGenome size: {genome_size:,} bp")
    print(f"Number of chunks: {len(chunks):,}")

    print("\nGenerating null distribution for p-value calculation...")
    background_seqs = []
    for record in SeqIO.parse(genome_file, "fasta"):
        seq = str(record.seq)
        for _ in range(100):
            start = np.random.randint(0, len(seq) - 20)
            background_seqs.append(seq[start:start+20])
    
    predictor = SinRBindingPredictor()
    predictor.generate_null_distribution(background_seqs)

    print("\nStep 5: Analyzing Binding Sites")
    print("-" * 30)
    print("Processing chunks in parallel...")
    process_func = partial(process_sequence_chunk, predictor=predictor)

    with Pool(processes=num_processes) as pool:
        results = list(tqdm(
            pool.imap(process_func, chunks),
            total=len(chunks),
            desc="Analyzing",
            ncols=100
        ))

    print("\nStep 6: Processing Results")
    print("-" * 30)
    flat_results = [item for sublist in results for item in sublist]
    final_results = pd.DataFrame(flat_results)
    final_results.drop_duplicates(subset=['position', 'sequence', 'contig'], inplace=True)

    if 'pvalue' in final_results.columns:
        pvals = final_results['pvalue'].values
        _, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')
        final_results['pvalue_corrected'] = pvals_corrected
        final_results['significant'] = final_results['pvalue_corrected'] <= 0.05

    final_results = final_results.sort_values('score', ascending=False)
    output_file = os.path.join(output_dir, f"{protein_name}_{genome_name}_binding_sites.csv")
    final_results.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")

    print("\nStep 7: Creating Visualizations")
    print("-" * 30)

    plt.figure(figsize=(10, 6))
    sns.histplot(data=final_results, x='score', bins=50)
    plt.title(f'{protein_name} Binding Score Distribution')
    plt.xlabel('Binding Score')
    plt.ylabel('Count')
    plot_file = os.path.join(output_dir, f"{protein_name}_{genome_name}_score_distribution.png")
    plt.savefig(plot_file)
    plt.close()

    plt.figure(figsize=(15, 6))
    if 'significant' in final_results.columns:
        sns.scatterplot(data=final_results, x='position', y='score',
                        hue='significant', palette={True: 'red', False: 'blue'},
                        alpha=0.5)
        plt.legend(title='Significant (FDR 5%)')
    else:
        sns.scatterplot(data=final_results, x='position', y='score', alpha=0.5)
    plt.title(f'{protein_name} Binding Sites Across Genome')
    plt.xlabel('Genome Position')
    plt.ylabel('Binding Score')
    position_plot_file = os.path.join(output_dir, f"{protein_name}_{genome_name}_genome_position.png")
    plt.savefig(position_plot_file)
    plt.close()

    print("\nStep 8: Summary")
    print("-" * 30)
    print(f"Total binding sites analyzed: {len(final_results):,}")
    if 'significant' in final_results.columns:
        print(f"Significant binding sites (FDR 5%): {final_results['significant'].sum():,}")
    print(f"\nTop 10 binding sites:")
    print(final_results.head(10).to_string())
    print(f"\nResults directory: {output_dir}")
    print("\nAnalysis complete!")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during analysis: {e}")
        sys.exit(1)
