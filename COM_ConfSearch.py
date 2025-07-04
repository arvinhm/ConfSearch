#!/usr/bin/env python3
"""
Optimized .com File Generator for Maximum Conformational Diversity
Creates ONE optimized .com file per peptide length (3-12 mers)
Parameters chosen to handle maximum conformational variation for AI training
Based on our earlier conversation about 2M diverse peptides for fragment-mimicry prediction
"""

import os
from pathlib import Path

def create_optimized_com_files():
    """
    Create ONE optimized .com file for each peptide length (3-12 mers)
    Parameters optimized for maximum conformational diversity and variation
    """
    
    print("🧬 Creating OPTIMIZED .com files for maximum conformational diversity")
    print("📚 For training AI on fragment-peptide mimicry prediction")
    print("=" * 70)
    
    # Create output directory
    output_dir = Path("OPTIMIZED_COM_FILES")
    output_dir.mkdir(exist_ok=True)
    
    # Optimized parameters for maximum conformational diversity per length
    # Based on literature review and our earlier discussion
    length_configs = {
        3: {
            'mcmm_steps': 200, 'mcop_operations': 2000, 'energy_window': 15.0, 'demx_steps': 15000,
            'method': 'mcmm', 'description': 'Short peptides - complete sampling possible'
        },
        4: {
            'mcmm_steps': 300, 'mcop_operations': 3000, 'energy_window': 18.0, 'demx_steps': 18000,
            'method': 'mcmm', 'description': 'Short peptides - enhanced sampling'
        },
        5: {
            'mcmm_steps': 500, 'mcop_operations': 5000, 'energy_window': 20.0, 'demx_steps': 20000,
            'method': 'mcmm', 'description': 'Medium peptides - increased diversity window'
        },
        6: {
            'mcmm_steps': 750, 'mcop_operations': 7500, 'energy_window': 22.0, 'demx_steps': 22000,
            'method': 'mcmm', 'description': 'Medium peptides - transition to complex sampling'
        },
        7: {
            'mcmm_steps': 1000, 'mcop_operations': 10000, 'energy_window': 25.0, 'demx_steps': 25000,
            'method': 'mixed', 'description': 'Long peptides - mixed MCMM+LMOD for diversity'
        },
        8: {
            'mcmm_steps': 1500, 'mcop_operations': 12500, 'energy_window': 28.0, 'demx_steps': 28000,
            'method': 'mixed', 'description': 'Long peptides - enhanced mixed sampling'
        },
        9: {
            'mcmm_steps': 2000, 'mcop_operations': 15000, 'energy_window': 30.0, 'demx_steps': 30000,
            'method': 'mixed', 'description': 'Very long peptides - comprehensive sampling'
        },
        10: {
            'mcmm_steps': 2500, 'mcop_operations': 17500, 'energy_window': 32.0, 'demx_steps': 32000,
            'method': 'mixed', 'description': 'Very long peptides - maximum diversity'
        },
        11: {
            'mcmm_steps': 3000, 'mcop_operations': 20000, 'energy_window': 35.0, 'demx_steps': 35000,
            'method': 'mixed_intensive', 'description': 'Extra long peptides - intensive mixed sampling'
        },
        12: {
            'mcmm_steps': 3500, 'mcop_operations': 22500, 'energy_window': 38.0, 'demx_steps': 38000,
            'method': 'mixed_intensive', 'description': 'Extra long peptides - maximum conformational coverage'
        }
    }
    
    created_files = []
    
    print("📋 Parameter Summary:")
    print("Length | Method           | MCMM | MCOP  | Energy | Runtime Est.")
    print("-------|------------------|------|-------|--------|------------")
    
    # Create optimized .com file for each length
    for length in range(3, 13):  # 3 to 12
        config = length_configs[length]
        
        # File names
        input_file = f"peptide_{length}mer.mae"
        output_file = f"peptide_{length}mer_conf_out.maegz"
        com_filename = f"peptide_{length}mer_optimized.com"
        com_filepath = output_dir / com_filename
        
        # Build .com content based on method
        com_content = f"""{input_file}
{output_file}
 DEBG       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 SEED       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 SOLV       3      1      0      0     0.0000     0.0000     0.0000     0.0000
 EXNB       0      5      0      0     8.0000    20.0000     4.0000     8.0000
 BDCO       0      0      0      0    89.4427 99999.0000     0.0000     0.0000
 FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000
 READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000"""

        # Method-specific parameters for optimal conformational diversity
        if config['method'] == 'mcmm':
            # Standard MCMM for shorter peptides
            com_content += f"""
 MCMM    {config['mcmm_steps']:4d}      0      0      0     0.0000     0.2500     0.0000     0.0000
 MCSS       2      0      0      0    50.0000     0.0000     0.0000     0.0000
 MCNV       1      0      0      0     0.0000     0.0000     0.0000     0.0000
 MCOP    {config['mcop_operations']:4d}      0      0      0     0.0000   {config['energy_window']:.1f}     0.0000     0.0000"""
        
        elif config['method'] == 'mixed':
            # Mixed MCMM + LMOD for longer peptides (better conformational sampling)
            mcmm_split = config['mcmm_steps'] // 2
            mcop_split = config['mcop_operations'] // 2
            com_content += f"""
 MCMM    {mcmm_split:4d}      0      0      0     0.0000     0.2500     0.0000     0.0000
 MCSS       2      0      0      0    50.0000     0.0000     0.0000     0.0000
 MCNV       1      0      0      0     0.0000     0.0000     0.0000     0.0000
 MCOP    {mcop_split:4d}      0      0      0     0.0000   {config['energy_window']:.1f}     0.0000     0.0000
 LMCS    {mcmm_split:4d}      0      0      0     0.0000     0.2500     0.0000     0.0000"""
        
        elif config['method'] == 'mixed_intensive':
            # Intensive mixed sampling for very long peptides
            mcmm_split = config['mcmm_steps'] // 3
            mcop_split = config['mcop_operations'] // 3
            com_content += f"""
 MCMM    {mcmm_split:4d}      0      0      0     0.0000     0.2500     0.0000     0.0000
 MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000
 MCNV       1      0      0      0     0.0000     0.0000     0.0000     0.0000
 MCOP    {mcop_split:4d}      0      0      0     0.0000   {config['energy_window']:.1f}     0.0000     0.0000
 LMCS    {mcmm_split:4d}      0      0      0     0.0000     0.2500     0.0000     0.0000
 LLMN    {mcmm_split:4d}      0      0      0     0.0000     0.2500     0.0000     0.0000"""

        # Common ending parameters optimized for diversity
        com_content += f"""
 DEMX       0  {config['demx_steps']:5d}      0      0   200.0000   100.0000     0.0000     0.0000
 MSYM       1      0      0      0     0.0000     0.0000     0.0000     0.0000
 AUOP       0      0      0      0    50.0000     0.0000     0.0000     0.0000
 AUTO       0      2      1      1     1.0000     1.0000     5.0000     1.0000
 CRMS       0      0      0      0     1.0000     0.0100    60.0000     0.0000
 NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000
 MINI       1      0   5000      0     0.0000     0.0010     0.0000     0.0000"""
        
        # Write the file
        with open(com_filepath, 'w') as f:
            f.write(com_content)
        
        created_files.append(com_filename)
        
        # Calculate runtime estimate
        runtime_est = config['mcop_operations'] / 1000 * 2  # rough estimate in hours
        
        print(f"  {length:2d}   | {config['method']:15s} | {config['mcmm_steps']:4d} | {config['mcop_operations']:5d} | {config['energy_window']:6.1f} | {runtime_est:4.0f}-{runtime_est*2:.0f}h")
    
    # Summary
    print("\n" + "=" * 70)
    print(f"🎉 SUCCESS! Created {len(created_files)} OPTIMIZED .com files")
    print("=" * 70)
    print("📁 Files created:")
    for filename in created_files:
        print(f"   ✅ {filename}")
    
    print(f"\n📂 Directory: {output_dir.absolute()}")
    print("\n🧬 OPTIMIZATION FEATURES:")
    print("   • Scaled parameters by peptide length")
    print("   • MCMM for short peptides (3-6 mers)")
    print("   • Mixed MCMM+LMOD for long peptides (7-10 mers)")  
    print("   • Intensive mixed sampling for very long peptides (11-12 mers)")
    print("   • Wide energy windows for conformational diversity")
    print("   • Optimized for your 2M peptide AI training dataset")
    
    print("\n🚀 READY FOR:")
    print("   • Maximum conformational diversity sampling")
    print("   • AI training on fragment-peptide mimicry")
    print("   • Handling diverse peptide sequences")
    print("   • Comprehensive conformational search")

if __name__ == "__main__":
    create_optimized_com_files()