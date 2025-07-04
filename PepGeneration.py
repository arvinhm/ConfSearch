#!/usr/bin/env python3
"""
Refined AI-Optimized Peptide Library Generator
Generates 2 million diverse peptides (3-12 mers) for conformational AI training
Optimized for fragment-peptide mimicry prediction with separate output files
"""

import numpy as np
import random
import itertools
import gc
import os
from pathlib import Path
from collections import defaultdict, Counter
from typing import List, Dict, Set, Tuple, Optional
import time
import json

class RefinedPeptideLibraryGenerator:
    def __init__(self):
        # Enhanced amino acid properties with more detailed characteristics
        self.aa_properties = {
            'A': {'hydrophobicity': 1.8, 'charge': 0, 'size': 67, 'flexibility': 0.7, 'beta_propensity': 0.9, 'alpha_propensity': 1.42, 'turn_propensity': 0.66},
            'C': {'hydrophobicity': 2.5, 'charge': 0, 'size': 86, 'flexibility': 0.4, 'beta_propensity': 1.1, 'alpha_propensity': 0.70, 'turn_propensity': 1.19},
            'D': {'hydrophobicity': -3.5, 'charge': -1, 'size': 91, 'flexibility': 0.8, 'beta_propensity': 0.7, 'alpha_propensity': 1.01, 'turn_propensity': 1.46},
            'E': {'hydrophobicity': -3.5, 'charge': -1, 'size': 109, 'flexibility': 0.8, 'beta_propensity': 0.7, 'alpha_propensity': 1.51, 'turn_propensity': 0.74},
            'F': {'hydrophobicity': 2.8, 'charge': 0, 'size': 135, 'flexibility': 0.5, 'beta_propensity': 1.2, 'alpha_propensity': 1.13, 'turn_propensity': 0.60},
            'G': {'hydrophobicity': -0.4, 'charge': 0, 'size': 48, 'flexibility': 1.0, 'beta_propensity': 0.6, 'alpha_propensity': 0.57, 'turn_propensity': 1.56},
            'H': {'hydrophobicity': -3.2, 'charge': 0.5, 'size': 118, 'flexibility': 0.6, 'beta_propensity': 0.9, 'alpha_propensity': 1.00, 'turn_propensity': 0.95},
            'I': {'hydrophobicity': 4.5, 'charge': 0, 'size': 124, 'flexibility': 0.3, 'beta_propensity': 1.4, 'alpha_propensity': 1.08, 'turn_propensity': 0.47},
            'K': {'hydrophobicity': -3.9, 'charge': 1, 'size': 135, 'flexibility': 0.9, 'beta_propensity': 0.8, 'alpha_propensity': 1.16, 'turn_propensity': 1.01},
            'L': {'hydrophobicity': 3.8, 'charge': 0, 'size': 124, 'flexibility': 0.4, 'beta_propensity': 1.2, 'alpha_propensity': 1.21, 'turn_propensity': 0.59},
            'M': {'hydrophobicity': 1.9, 'charge': 0, 'size': 124, 'flexibility': 0.5, 'beta_propensity': 1.0, 'alpha_propensity': 1.45, 'turn_propensity': 0.60},
            'N': {'hydrophobicity': -3.5, 'charge': 0, 'size': 96, 'flexibility': 0.7, 'beta_propensity': 0.7, 'alpha_propensity': 0.67, 'turn_propensity': 1.56},
            'P': {'hydrophobicity': -1.6, 'charge': 0, 'size': 90, 'flexibility': 0.1, 'beta_propensity': 0.5, 'alpha_propensity': 0.57, 'turn_propensity': 1.52},
            'Q': {'hydrophobicity': -3.5, 'charge': 0, 'size': 114, 'flexibility': 0.7, 'beta_propensity': 0.8, 'alpha_propensity': 1.11, 'turn_propensity': 0.98},
            'R': {'hydrophobicity': -4.5, 'charge': 1, 'size': 148, 'flexibility': 0.9, 'beta_propensity': 0.8, 'alpha_propensity': 0.98, 'turn_propensity': 0.95},
            'S': {'hydrophobicity': -0.8, 'charge': 0, 'size': 73, 'flexibility': 0.8, 'beta_propensity': 0.8, 'alpha_propensity': 0.77, 'turn_propensity': 1.43},
            'T': {'hydrophobicity': -0.7, 'charge': 0, 'size': 93, 'flexibility': 0.6, 'beta_propensity': 1.0, 'alpha_propensity': 0.83, 'turn_propensity': 0.96},
            'V': {'hydrophobicity': 4.2, 'charge': 0, 'size': 105, 'flexibility': 0.3, 'beta_propensity': 1.5, 'alpha_propensity': 1.06, 'turn_propensity': 0.50},
            'W': {'hydrophobicity': -0.9, 'charge': 0, 'size': 163, 'flexibility': 0.4, 'beta_propensity': 1.1, 'alpha_propensity': 1.08, 'turn_propensity': 0.96},
            'Y': {'hydrophobicity': -1.3, 'charge': 0, 'size': 141, 'flexibility': 0.5, 'beta_propensity': 1.2, 'alpha_propensity': 0.69, 'turn_propensity': 1.14}
        }
        
        self.amino_acids = list(self.aa_properties.keys())
        
        # Enhanced natural frequencies based on proteome analysis
        self.natural_freq = {
            'A': 0.083, 'C': 0.014, 'D': 0.054, 'E': 0.062, 'F': 0.039,
            'G': 0.068, 'H': 0.022, 'I': 0.056, 'K': 0.058, 'L': 0.096,
            'M': 0.024, 'N': 0.043, 'P': 0.051, 'Q': 0.040, 'R': 0.052,
            'S': 0.068, 'T': 0.058, 'V': 0.068, 'W': 0.013, 'Y': 0.032
        }
        
        # Enhanced functional motifs for comprehensive coverage
        self.functional_motifs = {
            'beta_turns': ['GP', 'GS', 'NG', 'DG', 'PG', 'DP', 'GG', 'PP', 'SG', 'GN'],
            'beta_strands': ['FY', 'VV', 'II', 'LV', 'FF', 'YY', 'WW', 'IV', 'VI', 'LI'],
            'charged_pairs': ['RK', 'DE', 'KE', 'RD', 'EK', 'DR', 'RR', 'DD', 'EE', 'KK'],
            'hydrophobic_clusters': ['LL', 'IV', 'FW', 'LI', 'VL', 'WF', 'FL', 'LF', 'VV', 'II'],
            'flexible_regions': ['GG', 'GA', 'AG', 'GN', 'NG', 'GS', 'SG', 'GD', 'AS', 'SA'],
            'aromatic_stacks': ['FF', 'FY', 'FW', 'YF', 'YY', 'YW', 'WF', 'WY', 'WW', 'HF'],
            'polar_networks': ['ST', 'SN', 'TN', 'SQ', 'TQ', 'NQ', 'SS', 'TT', 'NN', 'QQ'],
            'proline_hinges': ['PG', 'GP', 'PA', 'AP', 'PS', 'SP', 'PT', 'TP'],
            'disulfide_prone': ['CC', 'CG', 'GC', 'CS', 'SC', 'CT', 'TC'],
            'membrane_spans': ['LLL', 'III', 'VVV', 'FFF', 'LIV', 'VIL', 'FLI']
        }

        # Enhanced chemical property groups
        self.aa_groups = {
            'hydrophobic': ['I', 'L', 'V', 'F', 'W', 'Y', 'M', 'A'],
            'hydrophilic': ['R', 'K', 'D', 'E', 'N', 'Q', 'S', 'T', 'H'],
            'special': ['G', 'P', 'C'],
            'aromatic': ['F', 'W', 'Y', 'H'],
            'charged_pos': ['R', 'K', 'H'],
            'charged_neg': ['D', 'E'],
            'small': ['G', 'A', 'S', 'C'],
            'large': ['W', 'F', 'Y', 'R', 'K', 'Q', 'E'],
            'flexible': ['G', 'A', 'S', 'N', 'D', 'K', 'R'],
            'rigid': ['P', 'V', 'I', 'W'],
            'beta_forming': ['V', 'I', 'Y', 'F', 'L', 'W'],
            'alpha_forming': ['E', 'A', 'L', 'M', 'Q', 'K', 'R'],
            'turn_forming': ['G', 'N', 'P', 'S', 'D']
        }
        
        # Optimized distribution for 2 million peptides
        self.target_distribution = {
            3: 8000,        # All possible 3-mers (20^3 = 8,000)
            4: 160000,      # All possible 4-mers (20^4 = 160,000)
            5: 400000,      # Strategic 5-mers (from 3.2M possible)
            6: 350000,      # Strategic 6-mers
            7: 300000,      # Strategic 7-mers
            8: 250000,      # Strategic 8-mers
            9: 200000,      # Strategic 9-mers
            10: 150000,     # Strategic 10-mers
            11: 100000,     # Strategic 11-mers
            12: 82000       # Strategic 12-mers
        }
        # Total: 2,000,000 peptides

    def generate_all_short_peptides(self, length: int) -> List[str]:
        """Generate ALL possible peptides for lengths 3-4 with memory optimization"""
        print(f"Generating ALL {20**length:,} possible {length}-mers...")
        
        if length > 4:
            raise ValueError("Use sampling for peptides longer than 4 residues")
        
        # Memory-efficient generation using itertools
        all_peptides = []
        batch_size = 10000
        current_batch = []
        
        total_expected = 20 ** length
        processed = 0
        
        for combination in itertools.product(self.amino_acids, repeat=length):
            peptide = ''.join(combination)
            current_batch.append(peptide)
            processed += 1
            
            if len(current_batch) >= batch_size:
                all_peptides.extend(current_batch)
                current_batch = []
                
                if processed % 100000 == 0:
                    print(f"  Progress: {processed:,}/{total_expected:,} ({processed/total_expected*100:.1f}%)")
        
        # Add remaining peptides
        if current_batch:
            all_peptides.extend(current_batch)
        
        print(f"Generated {len(all_peptides):,} complete {length}-mers")
        return all_peptides

    def calculate_peptide_properties(self, sequence: str) -> Dict:
        """Calculate comprehensive chemical properties of a peptide"""
        props = {prop: [] for prop in ['hydrophobicity', 'charge', 'size', 'flexibility', 
                                     'beta_propensity', 'alpha_propensity', 'turn_propensity']}
        
        for aa in sequence:
            for prop in props:
                props[prop].append(self.aa_properties[aa][prop])
        
        # Calculate derived properties
        net_charge = sum(props['charge'])
        avg_hydrophobicity = np.mean(props['hydrophobicity'])
        charge_density = abs(net_charge) / len(sequence)
        
        return {
            'hydrophobicity': avg_hydrophobicity,
            'charge': net_charge,
            'charge_density': charge_density,
            'size': np.mean(props['size']),
            'flexibility': np.mean(props['flexibility']),
            'beta_propensity': np.mean(props['beta_propensity']),
            'alpha_propensity': np.mean(props['alpha_propensity']),
            'turn_propensity': np.mean(props['turn_propensity']),
            'length': len(sequence),
            'hydrophobic_moment': self._calculate_hydrophobic_moment(sequence),
            'amphipathicity': self._calculate_amphipathicity(sequence)
        }

    def _calculate_hydrophobic_moment(self, sequence: str) -> float:
        """Calculate hydrophobic moment for membrane interaction prediction"""
        if len(sequence) < 3:
            return 0.0
        
        angle_per_residue = 100  # degrees for alpha helix
        hydrophobic_sum_x = 0
        hydrophobic_sum_y = 0
        
        for i, aa in enumerate(sequence):
            angle = np.radians(i * angle_per_residue)
            hydrophobicity = self.aa_properties[aa]['hydrophobicity']
            hydrophobic_sum_x += hydrophobicity * np.cos(angle)
            hydrophobic_sum_y += hydrophobicity * np.sin(angle)
        
        return np.sqrt(hydrophobic_sum_x**2 + hydrophobic_sum_y**2) / len(sequence)

    def _calculate_amphipathicity(self, sequence: str) -> float:
        """Calculate amphipathicity score"""
        if len(sequence) < 2:
            return 0.0
        
        hydrophobic_count = sum(1 for aa in sequence if aa in self.aa_groups['hydrophobic'])
        hydrophilic_count = sum(1 for aa in sequence if aa in self.aa_groups['hydrophilic'])
        
        total = len(sequence)
        hydrophobic_fraction = hydrophobic_count / total
        hydrophilic_fraction = hydrophilic_count / total
        
        # Amphipathicity is maximized when both fractions are around 0.5
        return 4 * hydrophobic_fraction * hydrophilic_fraction

    def is_valid_sequence(self, sequence: str) -> bool:
        """Enhanced quality filters with stricter criteria for longer peptides"""
        if len(sequence) <= 4:
            return True  # Accept all short peptides
            
        # Count amino acid frequencies
        aa_counts = Counter(sequence)
        total_length = len(sequence)
        
        # Filter 1: No more than 35% of any single amino acid (stricter)
        for count in aa_counts.values():
            if count / total_length > 0.35:
                return False
        
        # Filter 2: No more than 25% proline (very strict)
        if aa_counts.get('P', 0) / total_length > 0.25:
            return False
            
        # Filter 3: Balanced hydrophobic content (20-70%)
        hydrophobic_count = sum(aa_counts.get(aa, 0) for aa in self.aa_groups['hydrophobic'])
        hydrophobic_fraction = hydrophobic_count / total_length
        if hydrophobic_fraction < 0.2 or hydrophobic_fraction > 0.7:
            return False
            
        # Filter 4: No long poly-X stretches (>3 consecutive for most, >2 for P,G)
        for i in range(len(sequence) - 2):
            if sequence[i] == sequence[i+1] == sequence[i+2]:
                if sequence[i] in ['P', 'G']:
                    return False  # No triple P or G
                if i < len(sequence) - 3 and sequence[i] == sequence[i+3]:
                    return False  # No quad repeats for any AA
        
        # Filter 5: Must have some structural diversity
        if len(set(sequence)) < max(2, len(sequence) // 4):
            return False
        
        # Filter 6: Reasonable charge distribution
        net_charge = sum(self.aa_properties[aa]['charge'] for aa in sequence)
        if abs(net_charge) > len(sequence) * 0.6:  # No more than 60% charged
            return False
            
        return True

    def generate_maximal_diversity_sequences(self, length: int, n_sequences: int) -> List[str]:
        """Advanced diversity generation using multiple strategies"""
        print(f"   Generating {n_sequences:,} diverse {length}-mers using maximal diversity approach...")
        
        sequences = set()
        strategies = [
            ('property_grid', 0.25),      # 25% - systematic property coverage
            ('motif_enhanced', 0.20),     # 20% - functional motifs
            ('chemical_extremes', 0.15),  # 15% - extreme properties  
            ('natural_sampling', 0.15),   # 15% - natural frequency
            ('fragment_based', 0.15),     # 15% - protein fragment mimics
            ('random_diverse', 0.10)      # 10% - pure random with diversity check
        ]
        
        for strategy, fraction in strategies:
            strategy_count = int(n_sequences * fraction)
            print(f"     {strategy}: targeting {strategy_count:,} peptides...")
            
            if strategy == 'property_grid':
                peptides = self._generate_property_grid_sequences(length, strategy_count)
            elif strategy == 'motif_enhanced':
                peptides = self._generate_motif_enhanced_sequences(length, strategy_count)
            elif strategy == 'chemical_extremes':
                peptides = self._generate_chemical_extremes(length, strategy_count)
            elif strategy == 'natural_sampling':
                peptides = self._generate_natural_frequency_sequences(length, strategy_count)
            elif strategy == 'fragment_based':
                peptides = self._generate_fragment_based_sequences(length, strategy_count)
            elif strategy == 'random_diverse':
                peptides = self._generate_random_diverse_sequences(length, strategy_count)
            
            # Add valid sequences
            valid_count = 0
            for peptide in peptides:
                if peptide and self.is_valid_sequence(peptide) and peptide not in sequences:
                    sequences.add(peptide)
                    valid_count += 1
                    if len(sequences) >= n_sequences:
                        break
            
            print(f"     {strategy}: added {valid_count:,} valid unique peptides")
            
            if len(sequences) >= n_sequences:
                break
        
        # Fill any remaining slots with high-quality random sequences
        while len(sequences) < n_sequences:
            peptide = self._generate_smart_random_sequence(length)
            if peptide and self.is_valid_sequence(peptide):
                sequences.add(peptide)
        
        result = list(sequences)[:n_sequences]
        print(f"   Completed: {len(result):,} unique diverse {length}-mers")
        return result

    def _generate_property_grid_sequences(self, length: int, n_sequences: int) -> List[str]:
        """Generate sequences covering systematic property grid"""
        sequences = []
        
        # Enhanced property grid with more dimensions
        hydrophob_points = np.linspace(-4.0, 4.5, 12)
        charge_points = list(range(-min(length//2, 6), min(length//2, 6) + 1))
        flex_points = np.linspace(0.2, 0.9, 8)
        amphipathicity_points = np.linspace(0.0, 1.0, 6)
        
        total_points = len(hydrophob_points) * len(charge_points) * len(flex_points) * len(amphipathicity_points)
        samples_per_point = max(1, n_sequences // total_points)
        
        for hydrophob in hydrophob_points:
            for charge in charge_points:
                for flex in flex_points:
                    for amphipathicity in amphipathicity_points:
                        if len(sequences) >= n_sequences:
                            break
                        
                        target_props = {
                            'hydrophobicity': hydrophob,
                            'charge': charge,
                            'flexibility': flex,
                            'amphipathicity': amphipathicity
                        }
                        
                        for _ in range(samples_per_point):
                            sequence = self._generate_property_targeted_sequence(length, target_props)
                            if sequence:
                                sequences.append(sequence)
                                if len(sequences) >= n_sequences:
                                    break
        
        return sequences[:n_sequences]

    def _generate_motif_enhanced_sequences(self, length: int, n_sequences: int) -> List[str]:
        """Generate sequences with embedded functional motifs"""
        sequences = []
        
        # Get suitable motifs for this length
        suitable_motifs = []
        for category, motifs in self.functional_motifs.items():
            for motif in motifs:
                if len(motif) <= length:
                    suitable_motifs.append((motif, category))
        
        if not suitable_motifs:
            return self._generate_natural_frequency_sequences(length, n_sequences)
        
        motif_usage = Counter()
        max_motif_usage = n_sequences // len(suitable_motifs) + 10
        
        attempts = 0
        while len(sequences) < n_sequences and attempts < n_sequences * 3:
            motif, category = random.choice(suitable_motifs)
            
            # Limit motif usage for diversity
            if motif_usage[motif] >= max_motif_usage:
                attempts += 1
                continue
            
            sequence = self._create_sequence_with_motif(motif, length)
            if sequence:
                sequences.append(sequence)
                motif_usage[motif] += 1
            
            attempts += 1
        
        return sequences

    def _generate_fragment_based_sequences(self, length: int, n_sequences: int) -> List[str]:
        """Generate sequences mimicking real protein fragments"""
        sequences = []
        
        # Common protein secondary structure patterns
        patterns = {
            'alpha_helix': [self.aa_groups['alpha_forming'], 0.6],
            'beta_strand': [self.aa_groups['beta_forming'], 0.5],
            'turn_region': [self.aa_groups['turn_forming'], 0.4],
            'flexible_loop': [self.aa_groups['flexible'], 0.5]
        }
        
        pattern_sequences = n_sequences // len(patterns)
        
        for pattern_name, (preferred_aas, preference_strength) in patterns.items():
            for _ in range(pattern_sequences):
                sequence = []
                for pos in range(length):
                    if random.random() < preference_strength:
                        aa = random.choice(preferred_aas)
                    else:
                        aa = random.choice(self.amino_acids)
                    sequence.append(aa)
                
                sequences.append(''.join(sequence))
                
                if len(sequences) >= n_sequences:
                    break
        
        return sequences

    def _generate_property_targeted_sequence(self, length: int, target_props: Dict) -> Optional[str]:
        """Generate sequence targeting specific properties with improved algorithm"""
        sequence = []
        
        # Calculate how many of each type of residue we need
        target_charge = target_props['charge']
        target_hydrophob = target_props['hydrophobicity']
        target_flex = target_props['flexibility']
        target_amphipathicity = target_props.get('amphipathicity', 0.5)
        
        # Phase 1: Place charged residues
        charge_positions = set()
        if target_charge > 0:
            charge_positions = set(random.sample(range(length), min(length, abs(target_charge))))
        elif target_charge < 0:
            charge_positions = set(random.sample(range(length), min(length, abs(target_charge))))
        
        # Phase 2: Determine hydrophobic/hydrophilic balance for amphipathicity
        remaining_positions = [i for i in range(length) if i not in charge_positions]
        hydrophobic_needed = int(len(remaining_positions) * target_amphipathicity)
        
        hydrophobic_positions = set(random.sample(remaining_positions, 
                                                min(len(remaining_positions), hydrophobic_needed)))
        
        # Phase 3: Build sequence
        for pos in range(length):
            weights = []
            for aa in self.amino_acids:
                props = self.aa_properties[aa]
                fitness = 1.0
                
                # Charge constraint
                if pos in charge_positions:
                    if target_charge > 0:
                        fitness *= 5.0 if props['charge'] > 0 else 0.1
                    else:
                        fitness *= 5.0 if props['charge'] < 0 else 0.1
                else:
                    fitness *= 2.0 if props['charge'] == 0 else 0.3
                
                # Hydrophobicity constraint
                if pos in hydrophobic_positions:
                    fitness *= 3.0 if aa in self.aa_groups['hydrophobic'] else 0.2
                else:
                    fitness *= 3.0 if aa in self.aa_groups['hydrophilic'] else 0.5
                
                # Flexibility preference
                flex_diff = abs(props['flexibility'] - target_flex)
                fitness *= 1.0 / (1.0 + flex_diff * 2)
                
                # Overall hydrophobicity preference
                hydrophob_diff = abs(props['hydrophobicity'] - target_hydrophob)
                fitness *= 1.0 / (1.0 + hydrophob_diff * 0.5)
                
                weights.append(max(fitness, 0.01))  # Minimum weight
            
            chosen_aa = random.choices(self.amino_acids, weights=weights, k=1)[0]
            sequence.append(chosen_aa)
        
        return ''.join(sequence)

    def _create_sequence_with_motif(self, motif: str, length: int) -> Optional[str]:
        """Create sequence with embedded motif"""
        if length == len(motif):
            return motif
        
        # Choose motif position
        max_start = length - len(motif)
        motif_start = random.randint(0, max_start)
        
        sequence = ['X'] * length
        
        # Insert motif
        for i, aa in enumerate(motif):
            sequence[motif_start + i] = aa
        
        # Fill remaining positions
        remaining_positions = [i for i, aa in enumerate(sequence) if aa == 'X']
        
        for pos in remaining_positions:
            # Use natural frequency with some bias towards sequence balance
            weights = [self.natural_freq[aa] for aa in self.amino_acids]
            
            # Slight bias for balance
            current_counts = Counter([aa for aa in sequence if aa != 'X'])
            for i, aa in enumerate(self.amino_acids):
                if current_counts[aa] > length * 0.3:  # Reduce over-represented AAs
                    weights[i] *= 0.5
            
            chosen_aa = random.choices(self.amino_acids, weights=weights, k=1)[0]
            sequence[pos] = chosen_aa
        
        return ''.join(sequence)

    def _generate_natural_frequency_sequences(self, length: int, n_sequences: int) -> List[str]:
        """Generate sequences using natural amino acid frequencies"""
        sequences = []
        weights = [self.natural_freq[aa] for aa in self.amino_acids]
        
        for _ in range(n_sequences):
            sequence = random.choices(self.amino_acids, weights=weights, k=length)
            sequences.append(''.join(sequence))
        
        return sequences

    def _generate_chemical_extremes(self, length: int, n_sequences: int) -> List[str]:
        """Generate sequences with extreme chemical properties"""
        sequences = []
        
        extreme_types = [
            ('max_hydrophobic', self.aa_groups['hydrophobic'], 0.7),
            ('max_hydrophilic', self.aa_groups['hydrophilic'], 0.7),
            ('max_positive', self.aa_groups['charged_pos'], 0.5),
            ('max_negative', self.aa_groups['charged_neg'], 0.5),
            ('max_aromatic', self.aa_groups['aromatic'], 0.4),
            ('max_small', self.aa_groups['small'], 0.6),
            ('max_large', self.aa_groups['large'], 0.6),
            ('max_flexible', self.aa_groups['flexible'], 0.7),
            ('max_rigid', self.aa_groups['rigid'], 0.5)
        ]
        
        sequences_per_type = n_sequences // len(extreme_types)
        
        for extreme_type, preferred_aas, fraction in extreme_types:
            for _ in range(sequences_per_type):
                sequence = []
                target_count = int(length * fraction)
                
                # Add preferred amino acids
                for _ in range(target_count):
                    sequence.append(random.choice(preferred_aas))
                
                # Fill remaining with diverse selection
                remaining = length - len(sequence)
                other_aas = [aa for aa in self.amino_acids if aa not in preferred_aas]
                for _ in range(remaining):
                    sequence.append(random.choice(other_aas))
                
                random.shuffle(sequence)
                sequences.append(''.join(sequence))
        
        return sequences

    def _generate_random_diverse_sequences(self, length: int, n_sequences: int) -> List[str]:
        """Generate random sequences with diversity checking"""
        sequences = []
        diversity_threshold = 0.3
        max_attempts = n_sequences * 5
        
        attempts = 0
        while len(sequences) < n_sequences and attempts < max_attempts:
            sequence = self._generate_smart_random_sequence(length)
            if sequence and self._is_diverse_enough(sequence, sequences, diversity_threshold):
                sequences.append(sequence)
            attempts += 1
        
        return sequences

    def _generate_smart_random_sequence(self, length: int) -> str:
        """Generate random sequence with smart constraints"""
        sequence = []
        
        # Ensure some diversity in the sequence
        min_unique = max(2, length // 3)
        used_aas = set()
        
        for i in range(length):
            if len(used_aas) < min_unique:
                # Force diversity
                available_aas = [aa for aa in self.amino_acids if aa not in used_aas]
                if available_aas:
                    aa = random.choice(available_aas)
                else:
                    aa = random.choice(self.amino_acids)
            else:
                # Normal selection with slight bias against over-representation
                current_counts = Counter(sequence)
                weights = []
                for aa in self.amino_acids:
                    weight = self.natural_freq[aa]
                    if current_counts[aa] > length * 0.3:
                        weight *= 0.3  # Reduce over-represented
                    weights.append(weight)
                
                aa = random.choices(self.amino_acids, weights=weights, k=1)[0]
            
            sequence.append(aa)
            used_aas.add(aa)
        
        return ''.join(sequence)

    def _is_diverse_enough(self, sequence: str, existing_sequences: List[str], threshold: float) -> bool:
        """Check if sequence is diverse enough compared to existing ones"""
        if not existing_sequences:
            return True
        
        # Sample existing sequences for efficiency
        sample_size = min(100, len(existing_sequences))
        sample = random.sample(existing_sequences, sample_size)
        
        for existing in sample:
            similarity = self._calculate_sequence_similarity(sequence, existing)
            if similarity > (1.0 - threshold):
                return False
        
        return True

    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity"""
        if len(seq1) != len(seq2):
            return 0.0
        
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1)

    def generate_comprehensive_library(self, target_total: int = 2_000_000, output_dir: str = "peptide_library") -> Dict:
        """Generate the complete comprehensive library with progress tracking"""
        print("🧬 COMPREHENSIVE PEPTIDE LIBRARY GENERATOR FOR AI TRAINING")
        print("=" * 80)
        print(f"Target: {target_total:,} diverse peptides (3-12 mers)")
        print("Strategy: Complete 3-4 mer coverage + maximal diversity 5-12 mers")
        print("Output: Separate files for each length")
        print("=" * 80)
        
        # Create output directory
        Path(output_dir).mkdir(exist_ok=True)
        
        all_statistics = {}
        generation_start_time = time.time()
        
        # PHASE 1: Complete short peptides
        print(f"\n🔬 PHASE 1: COMPLETE SHORT PEPTIDE COVERAGE")
        print("-" * 60)
        
        for length in [3, 4]:
            print(f"\nGenerating ALL {length}-mers...")
            start_time = time.time()
            
            peptides = self.generate_all_short_peptides(length)
            
            # Save to file immediately to manage memory
            filename = f"{output_dir}/peptides_length_{length}.txt"
            self._save_peptides_to_file(peptides, filename, length, "Complete enumeration")
            
            # Calculate statistics
            end_time = time.time()
            all_statistics[length] = {
                'count': len(peptides),
                'coverage': '100% (ALL POSSIBLE)',
                'strategy': 'Complete enumeration',
                'time_minutes': (end_time - start_time) / 60,
                'filename': filename
            }
            
            print(f"✅ {length}-mers: {len(peptides):,} peptides saved to {filename}")
            print(f"   Time: {(end_time - start_time)/60:.1f} minutes")
            
            # Clear memory
            del peptides
            gc.collect()
        
        # PHASE 2: Strategic longer peptides
        print(f"\n🎯 PHASE 2: STRATEGIC LONGER PEPTIDE SAMPLING")
        print("-" * 60)
        
        for length in range(5, 13):  # 5-12 mers
            target_count = self.target_distribution[length]
            print(f"\nGenerating {target_count:,} peptides of length {length}...")
            start_time = time.time()
            
            peptides = self.generate_maximal_diversity_sequences(length, target_count)
            
            # Save to file
            filename = f"{output_dir}/peptides_length_{length}.txt"
            self._save_peptides_to_file(peptides, filename, length, "Maximal diversity sampling")
            
            # Calculate statistics
            end_time = time.time()
            sample_props = self._calculate_sample_statistics(peptides[:1000])  # Sample for stats
            
            all_statistics[length] = {
                'count': len(peptides),
                'coverage': f'{len(peptides)/(20**length)*100:.8f}%',
                'strategy': 'Maximal diversity sampling',
                'time_minutes': (end_time - start_time) / 60,
                'filename': filename,
                **sample_props
            }
            
            print(f"✅ {length}-mers: {len(peptides):,} peptides saved to {filename}")
            print(f"   Time: {(end_time - start_time)/60:.1f} minutes")
            print(f"   Coverage: {len(peptides)/(20**length)*100:.8f}% of theoretical space")
            
            # Clear memory
            del peptides
            gc.collect()
        
        # Generate comprehensive summary
        total_time = time.time() - generation_start_time
        total_peptides = sum(stats['count'] for stats in all_statistics.values())
        
        print(f"\n" + "=" * 80)
        print(f"🎉 LIBRARY GENERATION COMPLETE!")
        print(f"=" * 80)
        print(f"Total peptides: {total_peptides:,}")
        print(f"Target was: {target_total:,}")
        print(f"Achievement: {total_peptides/target_total*100:.1f}%")
        print(f"Total time: {total_time/60:.1f} minutes ({total_time/3600:.1f} hours)")
        print(f"Files saved in: {output_dir}/")
        
        # Save master summary
        self._save_master_summary(all_statistics, output_dir, total_peptides, total_time)
        
        return all_statistics

    def _save_peptides_to_file(self, peptides: List[str], filename: str, length: int, strategy: str):
        """Save peptides to file with metadata header"""
        with open(filename, 'w') as f:
            f.write(f"# PEPTIDE LIBRARY - LENGTH {length}\n")
            f.write(f"# ================================\n")
            f.write(f"# Count: {len(peptides):,} peptides\n")
            f.write(f"# Strategy: {strategy}\n")
            f.write(f"# Generated for AI conformational training\n")
            f.write(f"# Format: One peptide sequence per line\n")
            f.write(f"#\n")
            
            for peptide in sorted(peptides):  # Sort for reproducibility
                f.write(f"{peptide}\n")

    def _calculate_sample_statistics(self, sample_peptides: List[str]) -> Dict:
        """Calculate statistics for a sample of peptides"""
        if not sample_peptides:
            return {}
        
        properties = [self.calculate_peptide_properties(seq) for seq in sample_peptides]
        
        return {
            'avg_hydrophobicity': np.mean([p['hydrophobicity'] for p in properties]),
            'avg_charge': np.mean([p['charge'] for p in properties]),
            'avg_flexibility': np.mean([p['flexibility'] for p in properties]),
            'charge_range': (min([p['charge'] for p in properties]), 
                           max([p['charge'] for p in properties])),
            'hydrophob_range': (min([p['hydrophobicity'] for p in properties]), 
                              max([p['hydrophobicity'] for p in properties])),
            'avg_amphipathicity': np.mean([p['amphipathicity'] for p in properties])
        }

    def _save_master_summary(self, all_statistics: Dict, output_dir: str, total_peptides: int, total_time: float):
        """Save comprehensive summary and statistics"""
        summary_file = f"{output_dir}/MASTER_SUMMARY.txt"
        
        with open(summary_file, 'w') as f:
            f.write("COMPREHENSIVE PEPTIDE LIBRARY - MASTER SUMMARY\n")
            f.write("=" * 60 + "\n")
            f.write(f"Generated for AI conformational training and fragment-peptide mimicry prediction\n\n")
            
            f.write(f"OVERVIEW:\n")
            f.write(f"Total peptides: {total_peptides:,}\n")
            f.write(f"Generation time: {total_time/60:.1f} minutes ({total_time/3600:.1f} hours)\n")
            f.write(f"Peptide lengths: 3-12 amino acids\n")
            f.write(f"Files generated: {len(all_statistics)} separate files\n\n")
            
            f.write(f"LENGTH DISTRIBUTION:\n")
            f.write("-" * 40 + "\n")
            
            for length in sorted(all_statistics.keys()):
                stats = all_statistics[length]
                f.write(f"Length {length:2d}: {stats['count']:,} peptides ")
                f.write(f"({stats['coverage']}) - {stats['filename']}\n")
                f.write(f"           Strategy: {stats['strategy']}\n")
                f.write(f"           Time: {stats['time_minutes']:.1f} min\n")
                
                if 'avg_hydrophobicity' in stats:
                    f.write(f"           Avg hydrophobicity: {stats['avg_hydrophobicity']:.2f}\n")
                    f.write(f"           Avg charge: {stats['avg_charge']:.1f}\n")
                    f.write(f"           Charge range: {stats['charge_range']}\n")
                    f.write(f"           Hydrophobicity range: {stats['hydrophob_range']}\n")
                
                f.write("\n")
            
            f.write(f"USAGE INSTRUCTIONS:\n")
            f.write("-" * 40 + "\n")
            f.write("1. Each file contains peptides of a single length\n")
            f.write("2. Files are ready for MacroModel conformational search\n")
            f.write("3. Use appropriate .com file parameters for each length\n")
            f.write("4. Combine results for comprehensive AI training dataset\n")
        
        # Also save as JSON for programmatic access
        json_file = f"{output_dir}/library_statistics.json"
        with open(json_file, 'w') as f:
            json.dump(all_statistics, f, indent=2, default=str)
        
        print(f"📊 Master summary saved to: {summary_file}")
        print(f"📈 Statistics JSON saved to: {json_file}")

# Main execution function
def main():
    """Main execution with user-friendly interface"""
    print("🚀 Starting Comprehensive Peptide Library Generation...")
    
    generator = RefinedPeptideLibraryGenerator()
    
    # Generate library
    statistics = generator.generate_comprehensive_library(
        target_total=2_000_000,
        output_dir="comprehensive_peptide_library"
    )
    
    print(f"\n✨ SUCCESS! Comprehensive peptide library ready for AI training!")
    print(f"📁 Check the 'comprehensive_peptide_library' directory for all files")
    print(f"🧬 {sum(stats['count'] for stats in statistics.values()):,} unique peptides generated")
    print(f"💡 Each length has its own file for easy MacroModel processing")

if __name__ == "__main__":
    main()