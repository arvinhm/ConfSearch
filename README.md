# Peptide Conformational Search Library

🧬 **Automated MacroModel Conformational Search Pipeline for Peptide Diversity Generation**

A comprehensive Python toolkit for generating maximum conformational diversity in peptides (3-12 amino acids) using Schrödinger's MacroModel. Optimized for AI training datasets and fragment-peptide mimicry prediction research.

## 🎯 Purpose

This project generates diverse peptide conformations for training AI models on fragment-peptide mimicry prediction. The pipeline creates optimized MacroModel configurations tailored to each peptide length, ensuring maximum conformational sampling for comprehensive dataset generation.

## ✨ Key Features

- **🔧 Intelligent Parameter Scaling**: Automatically adjusts MacroModel parameters based on peptide length (3-12 mers)
- **⚡ Multi-Core Acceleration**: Leverages all CPU cores for faster conformational searches
- **📊 Batch Processing**: Run multiple conformational searches concurrently
- **🎛️ Method Optimization**: 
  - MCMM sampling for short peptides (3-6 mers)
  - Mixed MCMM+LMOD for medium peptides (7-10 mers)
  - Intensive mixed sampling for long peptides (11-12 mers)
- **🤖 AI-Ready Output**: Optimized for training fragment-mimicry prediction models
- **📁 Organized Structure**: Systematic organization by peptide length with comprehensive libraries

## 🏗️ Project Structure

```
Peptide_ConfSearch/
├── COM_ConfSearch.py              # Generate optimized .com files
├── ConfSearch.py                  # Automated MacroModel runner
├── OPTIMIZED_COM_FILES/           # Generated .com configurations
│   ├── peptide_3mer_optimized.com
│   ├── peptide_4mer_optimized.com
│   └── ... (up to 12mer)
├── comprehensive_peptide_library/  # Complete peptide library
│   ├── 3AA/                       # 3 amino acid peptides
│   │   ├── build_peptide.py
│   │   ├── confgen.py
│   │   ├── ConfSearch.py
│   │   ├── sequences.txt
│   │   └── ...
│   ├── 4AA/ ... 12AA/             # 4-12 amino acid peptides
│   └── Text_Lib/                  # Peptide sequence libraries
│       ├── library_statistics.json
│       ├── MASTER_SUMMARY.txt
│       └── peptides_length_*.txt
```

## 🛠️ Requirements

### Software Dependencies
- **Schrödinger Suite** (2024.1 or later)
  - MacroModel for conformational searches
  - Maestro for structure preparation
- **Python 3.7+**
- **Windows 10/11** (optimized for Windows installations)

### Python Packages
```bash
# Standard library only - no additional packages required
```

## 📥 Installation

1. **Install Schrödinger Suite**
   - Download from [Schrödinger website](https://www.schrodinger.com/)
   - Install in standard location (e.g., `C:\Program Files\Schrodinger2025-2`)

2. **Set Environment Variable** (Optional)
   ```powershell
   # Windows PowerShell
   $env:SCHRODINGER = "C:\Program Files\Schrodinger2025-2"
   ```

3. **Clone Repository**
   ```bash
   git clone https://github.com/yourusername/Peptide_ConfSearch.git
   cd Peptide_ConfSearch
   ```

## 🚀 Quick Start

### 1. Generate Optimized Configuration Files

```bash
python COM_ConfSearch.py
```

This creates optimized `.com` files for all peptide lengths (3-12 mers) with parameters tuned for maximum conformational diversity.

**Output:**
```
🧬 Creating OPTIMIZED .com files for maximum conformational diversity
📚 For training AI on fragment-peptide mimicry prediction
======================================================================

📋 Parameter Summary:
Length | Method           | MCMM | MCOP  | Energy | Runtime Est.
-------|------------------|------|-------|--------|------------
   3   | mcmm             |  200 |  2000 |   15.0 |    4-8h
   4   | mcmm             |  300 |  3000 |   18.0 |    6-12h
   5   | mcmm             |  500 |  5000 |   20.0 |   10-20h
  ...  | ...              | ...  |  ...  |   ...  |   ...
  12   | mixed_intensive  | 3500 | 22500 |   38.0 |   45-90h
```

### 2. Run Conformational Searches

#### Interactive Mode
```bash
python ConfSearch.py
```

#### Single Job
```bash
python ConfSearch.py peptide_6mer_optimized.com
```

#### Batch Processing
```bash
python ConfSearch.py
# Select option 2 for batch processing
```

## 📖 Detailed Usage

### Configuration Generation (`COM_ConfSearch.py`)

Creates length-specific `.com` files with optimized parameters:

- **Energy Windows**: Scaled from 15.0 kcal/mol (3-mers) to 38.0 kcal/mol (12-mers)
- **MCMM Steps**: Scaled from 200 (3-mers) to 3,500 (12-mers)
- **Sampling Methods**:
  - **Short peptides** (3-6): Pure MCMM for complete sampling
  - **Medium peptides** (7-10): Mixed MCMM+LMOD for enhanced diversity
  - **Long peptides** (11-12): Intensive mixed sampling for maximum coverage

### Conformational Search Runner (`ConfSearch.py`)

#### Features:
- **Auto-detection** of Schrödinger installation
- **Multi-core support** with automatic CPU detection
- **Real-time output** monitoring
- **Batch processing** with concurrent job management
- **Error handling** and progress tracking

#### Command Line Options:
```bash
# Run with specific core count
python ConfSearch.py peptide_8mer_optimized.com 8

# Interactive mode with core selection
python ConfSearch.py
```

## 🧪 Example Workflow

### Complete Pipeline for 6-mer Peptides

1. **Prepare Input Structure**
   ```bash
   cd comprehensive_peptide_library/6AA/
   # Ensure peptide_6mer.mae exists
   ```

2. **Generate Optimized Configuration**
   ```bash
   python ../../COM_ConfSearch.py
   # Creates peptide_6mer_optimized.com
   ```

3. **Run Conformational Search**
   ```bash
   python ../../ConfSearch.py peptide_6mer_optimized.com
   ```

4. **Expected Output**
   ```
   📁 Output files created:
      ✅ peptide_6mer_conf_out.maegz (150.32 MB)
      ✅ peptide_6mer_optimized.log (2.45 MB)
   ```

## 🎛️ Configuration Details

### Parameter Optimization Strategy

| Peptide Length | Primary Method | MCMM Steps | Energy Window | Expected Runtime |
|---------------|----------------|------------|---------------|------------------|
| 3-6 mers      | MCMM          | 200-750    | 15.0-22.0     | 4-15 hours       |
| 7-10 mers     | Mixed         | 1000-2500  | 25.0-32.0     | 20-35 hours      |
| 11-12 mers    | Intensive     | 3000-3500  | 35.0-38.0     | 60-90 hours      |

### Energy Window Rationale
- **Wide energy windows** ensure capture of high-energy but biologically relevant conformations
- **Gradual scaling** accommodates increased conformational complexity with length
- **AI training optimization** prioritizes diversity over energy minimization

## 📊 Performance Expectations

### Computational Requirements
- **CPU**: Multi-core processor recommended (8+ cores optimal)
- **Memory**: 8GB+ RAM for longer peptides
- **Storage**: ~100-500 MB per peptide conformational search
- **Time**: 4-90 hours per peptide (length-dependent)

### Output Statistics
- **3-6 mers**: 50-200 unique conformations typically generated
- **7-10 mers**: 100-500 unique conformations typically generated
- **11-12 mers**: 200-1000+ unique conformations typically generated

## 🔬 Research Applications

### Primary Use Cases
- **Fragment-based drug design** conformational libraries
- **AI model training** for peptide-fragment mimicry prediction
- **Pharmacophore modeling** with diverse peptide conformations
- **Structure-activity relationship** studies

### AI Training Optimization
- Parameters specifically tuned for **maximum conformational diversity**
- **Wide energy windows** capture biologically relevant high-energy states
- **Systematic scaling** ensures consistent quality across peptide lengths
- **Batch processing** enables large-scale dataset generation

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Commit changes: `git commit -am 'Add feature'`
4. Push to branch: `git push origin feature-name`
5. Submit a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Schrödinger LLC** for MacroModel computational engine
- **Research community** for conformational search methodology development
- **Fragment-based drug design** field for motivating applications

## 📞 Support

For questions, issues, or contributions:
- 🐛 [Report bugs](https://github.com/yourusername/Peptide_ConfSearch/issues)
- 💡 [Request features](https://github.com/yourusername/Peptide_ConfSearch/issues)
- 📚 [View documentation](https://github.com/yourusername/Peptide_ConfSearch/wiki)

---

**Built for advancing AI-driven drug discovery through comprehensive peptide conformational analysis** 🧬🤖 
