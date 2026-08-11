"""
AI-Powered Genomic Insights
=============================

Generates expert-level DNA analysis insights using a local LLM
via Ollama. Provides graceful fallback when Ollama is not available.
"""

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)


def generate_ai_insights(
    sequence_name: str,
    stats: Dict,
    sequence: str,
    model: str = 'llama3.2',
) -> str:
    """
    Generate AI-powered insights about a DNA sequence using a local LLM.

    Args:
        sequence_name: Name of the sequence.
        stats: Dictionary from calculate_dna_stats.
        sequence: The actual DNA sequence string.
        model: Ollama model name (default: 'llama3.2').

    Returns:
        String with AI-generated insights, or a fallback analysis
        if Ollama is not available.
    """
    logger.info("Generating AI insights for '%s' using model '%s'...", sequence_name, model)

    prompt = f"""You are a genomics expert analyzing ancient DNA sequences. Provide a brief, professional analysis.

DNA Sequence Analysis:
- Sequence Name: {sequence_name}
- Length: {stats['length']} base pairs
- GC Content: {stats['gc_content']}%
- AT/GC Ratio: {stats.get('at_gc_ratio', 'N/A')}
- Base Composition:
  * Adenine (A): {stats['percentages']['A']}%
  * Thymine (T): {stats['percentages']['T']}%
  * Cytosine (C): {stats['percentages']['C']}%
  * Guanine (G): {stats['percentages']['G']}%

Sequence preview (first 100 bp): {sequence[:100]}

Provide a concise analysis covering:
1. Quality assessment based on GC content
2. Notable patterns or characteristics
3. Potential organism type hints (prokaryote vs eukaryote based on GC content)
4. Signs of ancient DNA damage (C→T deamination patterns)
5. Recommendations for further analysis

Keep response under 250 words and professional."""

    try:
        import ollama
        response = ollama.generate(model=model, prompt=prompt)
        return response['response']
    except ImportError:
        logger.warning("Ollama package not installed. Using fallback analysis.")
        return _fallback_analysis(sequence_name, stats)
    except Exception as e:
        logger.warning("Ollama unavailable: %s. Using fallback analysis.", e)
        return _fallback_analysis(sequence_name, stats)


def _fallback_analysis(sequence_name: str, stats: Dict) -> str:
    """Generate a rule-based fallback analysis when AI is not available."""
    gc = stats['gc_content']
    length = stats['length']

    # Organism type hint
    if gc < 30:
        organism_hint = "AT-rich — consistent with certain parasitic organisms or AT-rich genomes"
    elif gc < 45:
        organism_hint = "moderate GC — consistent with many eukaryotic organisms including mammals"
    elif gc < 55:
        organism_hint = "balanced GC — could be eukaryotic or prokaryotic"
    elif gc < 65:
        organism_hint = "GC-rich — common in certain prokaryotes (Actinobacteria)"
    else:
        organism_hint = "very high GC — unusual, may indicate thermophilic organisms"

    # Quality assessment
    if 35 <= gc <= 65:
        quality = "GC content is within normal range, suggesting good sequence quality."
    else:
        quality = "GC content is outside typical range — may indicate bias, contamination, or specialized genome."

    return f"""📊 Automated Analysis Report (AI Unavailable)

**Sequence:** {sequence_name}
**Length:** {length} bp | **GC Content:** {gc}%

**Quality Assessment:** {quality}

**Organism Hint:** {organism_hint}

**Base Balance:** AT/GC ratio is {stats.get('at_gc_ratio', 'N/A')}, purine/pyrimidine ratio is {stats.get('purine_pyrimidine_ratio', 'N/A')}.

**Recommendations:**
- Run sequence alignment against NCBI BLAST for species identification
- Check for C→T damage patterns characteristic of ancient DNA
- Perform contamination screening against modern reference databases

⚠️ *For AI-powered insights, install and start Ollama: `ollama serve` then `ollama pull llama3.2`*"""
