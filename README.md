# Ancient DNA Analyzer

A Python project for exploring DNA sequences, including sample data used for paleogenomics demonstrations. It accepts FASTA files or raw sequences and provides sequence statistics and basic analysis through a CLI or Streamlit dashboard. It does not authenticate a sample's age or origin, align reads, call variants from sequencing data, or validate historical or medical conclusions.

## What it does
- Parses FASTA files and validates DNA sequences; computes base composition, GC content and dinucleotide statistics.
- Explores codon use, open reading frames, protein translation, GC skew and restriction sites.
- Compares two sequences with a direct, position-by-position mutation summary. This is not a sequence aligner.
- Displays plots in Streamlit. Optional local Ollama output offers a narrative summary; without Ollama, it uses a rule-based fallback. Treat those narratives as exploratory, not scientific findings.

## Run locally
Requires Python 3.9+.

```bash
git clone https://github.com/rudra0812/ancient-dna-analyzer.git
cd ancient-dna-analyzer
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

For the CLI:

```bash
python cli.py analyze --fasta data/samples/sample_sequences.fasta --report
python cli.py analyze --sequence ATGCGATCGATCGATCG
python cli.py compare file1.fasta file2.fasta
```

Run tests with `python -m pytest tests/`. See `src/ancient_dna/` for the core, genomics, visualization, AI and reporting modules. Licensed under MIT; see `LICENSE`.
