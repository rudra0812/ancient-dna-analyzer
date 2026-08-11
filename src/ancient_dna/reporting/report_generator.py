"""
DNA Analysis Report Generator
===============================

Generates comprehensive analysis reports combining statistics,
visualizations, genomic engineering analysis, and AI insights
into a single pipeline.
"""

import logging
import os
from typing import Optional

from src.ancient_dna.core.validator import validate_dna_sequence
from src.ancient_dna.core.stats import calculate_dna_stats
from src.ancient_dna.ai.insights import generate_ai_insights

logger = logging.getLogger(__name__)


def generate_dna_report(
    sequence_name: str,
    sequence: str,
    save_report: bool = False,
    output_dir: str = ".",
    ai_model: str = "llama3.2",
) -> Optional[dict]:
    """
    Generate a complete analysis report with statistics and AI insights.

    Args:
        sequence_name: Name of the sequence.
        sequence: DNA sequence string.
        save_report: If True, saves report to a text file.
        output_dir: Directory to save reports in.
        ai_model: Ollama model name for AI insights.

    Returns:
        Dict with all analysis results, or None if validation fails.
    """
    separator = "=" * 70

    print(f"\n{separator}")
    print(f"🧬 DNA ANALYSIS REPORT: {sequence_name}")
    print(f"{separator}\n")

    # Step 1: Validate
    print("Step 1: Validating sequence...")
    try:
        validate_dna_sequence(sequence)
        print("  ✅ Sequence is valid")
    except ValueError as e:
        print(f"  ❌ Validation failed: {e}")
        return None

    # Step 2: Calculate statistics
    print("\nStep 2: Calculating statistics...")
    stats = calculate_dna_stats(sequence)
    if not stats:
        print("  ❌ Statistics calculation failed!")
        return None

    print(f"  ✅ Length: {stats['length']} bp")
    print(f"  ✅ GC Content: {stats['gc_content']}%")
    print(f"  ✅ AT/GC Ratio: {stats['at_gc_ratio']}")
    print(f"  ✅ Base counts: A={stats['counts']['A']}, T={stats['counts']['T']}, "
          f"C={stats['counts']['C']}, G={stats['counts']['G']}")

    # Step 3: Generate AI insights
    print("\nStep 3: Generating AI-powered insights...")
    ai_insights = generate_ai_insights(sequence_name, stats, sequence, model=ai_model)

    print(f"\n{separator}")
    print("🤖 AI-GENERATED INSIGHTS")
    print(separator)
    print(ai_insights)
    print(separator)

    # Step 4: Save report
    if save_report:
        os.makedirs(output_dir, exist_ok=True)
        report_filename = os.path.join(output_dir, f"{sequence_name.replace(' ', '_')}_report.txt")
        with open(report_filename, 'w', encoding='utf-8') as f:
            f.write(f"DNA ANALYSIS REPORT\n")
            f.write(f"{separator}\n\n")
            f.write(f"Sequence Name: {sequence_name}\n")
            f.write(f"Length: {stats['length']} bp\n")
            f.write(f"GC Content: {stats['gc_content']}%\n")
            f.write(f"AT/GC Ratio: {stats['at_gc_ratio']}\n")
            f.write(f"Purine/Pyrimidine Ratio: {stats['purine_pyrimidine_ratio']}\n\n")
            f.write(f"Base Composition:\n")
            for base, pct in stats['percentages'].items():
                f.write(f"  {base}: {pct}% ({stats['counts'][base]} bases)\n")
            f.write(f"\n{separator}\n")
            f.write(f"AI-GENERATED INSIGHTS\n")
            f.write(f"{separator}\n\n")
            f.write(ai_insights)

        full_path = os.path.abspath(report_filename)
        print(f"\n✅ Report saved: {full_path}")

    print(f"\n{separator}")
    print("✅ ANALYSIS COMPLETE!")
    print(f"{separator}\n")

    return {
        'sequence_name': sequence_name,
        'stats': stats,
        'ai_insights': ai_insights,
    }
