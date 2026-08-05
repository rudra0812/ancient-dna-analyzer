import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Dict, List
import os
import ollama

# Import the logic functions from our new module
from dna_logic import (calculate_dna_stats, read_fasta_file, transcribe_dna_to_rna,
                       validate_dna_sequence)

# Module-level constants
PLOT_COLORS: List[str] = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A']

# =============================================================================
# VISUALIZATION
# =============================================================================
def visualize_dna_stats(sequence_name: str, stats: Dict, save_plot: bool = False) -> None:
    """
    Create an interactive visualization of DNA statistics using Plotly.

    Args:
        sequence_name: Name of the sequence.
        stats: Dictionary from calculate_dna_stats.
        save_plot: If True, saves the interactive plot as an HTML file.
    """
    if not stats:
        print("No statistics to visualize.")
        return

    bases = list(stats["counts"].keys())
    counts = list(stats["counts"].values())
    percentages = list(stats["percentages"].values())

    # Create a figure with 1 row and 2 columns
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{'type': 'bar'}, {'type': 'pie'}]],
        subplot_titles=("Base Composition", f"Base Distribution (GC: {stats['gc_content']}%)")
    )

    # Subplot 1: Bar Chart
    fig.add_trace(go.Bar(
        x=bases,
        y=counts,
        marker_color=PLOT_COLORS,
        text=counts,
        textposition='auto',
        name='Counts'
    ), row=1, col=1)

    # Subplot 2: Pie Chart
    fig.add_trace(go.Pie(
        labels=bases,
        values=percentages,
        marker_colors=PLOT_COLORS,
        hoverinfo='label+percent',
        textinfo='percent+label'
    ), row=1, col=2)

    # Update layout for a polished look
    fig.update_layout(title_text=f"Analysis for: {sequence_name}", showlegend=False)
    fig.update_yaxes(title_text="Count", row=1, col=1)

    if save_plot:
        filename = f"{sequence_name.replace(' ', '_')}.html"
        fig.write_html(filename)
        print(f"Saved interactive plot: {filename}")

    fig.show()

# =============================================================================
# AI ANALYSIS & REPORTING
# =============================================================================
# Function 5: AI-Powered DNA Insights
def generate_ai_insights(sequence_name, stats, sequence):
    """
    Generate AI-powered insights about a DNA sequence using local LLM.
    
    Args:
        sequence_name: Name of the sequence
        stats: Dictionary from calculate_dna_stats
        sequence: The actual DNA sequence string
        
    Returns:
        String with AI-generated insights
    """
    print(f"\n🤖 Generating AI insights for {sequence_name}...")
    
    # Create a detailed prompt for the AI
    prompt = f"""You are a genomics expert analyzing DNA sequences. Provide a brief, professional analysis.

DNA Sequence Analysis:
- Sequence Name: {sequence_name}
- Length: {stats['length']} base pairs
- GC Content: {stats['gc_content']}%
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
4. Any recommendations for further analysis

Keep response under 200 words and professional."""

    try:
        # Call Ollama AI
        response = ollama.generate(
            model='llama3.2',
            prompt=prompt
        )
        
        return response['response']
        
    except Exception as e:
        return f"Error generating AI insights: {e}\nMake sure Ollama is running and llama3.2 is installed."


# Function 6: Generate Complete DNA Report
def generate_dna_report(sequence_name, sequence, save_report=False):
    """
    Generate a complete analysis report with statistics, visualization, and AI insights.
    
    Args:
        sequence_name: Name of the sequence
        sequence: DNA sequence string
        save_report: If True, saves report to text file
    """
    print(f"\n{'='*70}")
    print(f"DNA ANALYSIS REPORT: {sequence_name}")
    print(f"{'='*70}\n")
    
    # Step 1: Validate
    print("Step 1: Validating sequence...")
    if not validate_dna_sequence(sequence):
        print("❌ Sequence validation failed!")
        return
    
    # Step 2: Calculate statistics
    print("\nStep 2: Calculating statistics...")
    stats = calculate_dna_stats(sequence)
    
    if not stats:
        print("❌ Statistics calculation failed!")
        return
    
    print(f"✅ Length: {stats['length']} bp")
    print(f"✅ GC Content: {stats['gc_content']}%")
    print(f"✅ Base counts: A={stats['counts']['A']}, T={stats['counts']['T']}, "
          f"C={stats['counts']['C']}, G={stats['counts']['G']}")
    
    # Step 3: Generate AI insights
    print("\nStep 3: Generating AI-powered insights...")
    ai_insights = generate_ai_insights(sequence_name, stats, sequence)
    
    print("\n" + "="*70)
    print("🤖 AI-GENERATED INSIGHTS")
    print("="*70)
    print(ai_insights)
    print("="*70)
    
    # Step 4: Create visualization
    print("\nStep 4: Creating visualizations...")
    visualize_dna_stats(sequence_name, stats)
    
    # Step 5: Save report if requested
    if save_report:
        report_filename = f"{sequence_name.replace(' ', '_')}_report.txt"
        with open(report_filename, 'w') as f:
            f.write(f"DNA ANALYSIS REPORT\n")
            f.write(f"{'='*70}\n\n")
            f.write(f"Sequence Name: {sequence_name}\n")
            f.write(f"Length: {stats['length']} bp\n")
            f.write(f"GC Content: {stats['gc_content']}%\n\n")
            f.write(f"Base Composition:\n")
            for base, pct in stats['percentages'].items():
                f.write(f"  {base}: {pct}% ({stats['counts'][base]} bases)\n")
            f.write(f"\n{'='*70}\n")
            f.write(f"AI-GENERATED INSIGHTS\n")
            f.write(f"{'='*70}\n\n")
            f.write(ai_insights)
        
        full_path = os.path.abspath(report_filename)
        print(f"✅ Report saved as: {full_path}")
    
    print(f"\n{'='*70}")
    print("✅ ANALYSIS COMPLETE!")
    print(f"{'='*70}\n")

# ===============================
# MAIN EXECUTION BLOCK
# ===============================
if __name__ == "__main__":
    print("🧬 ANCIENT DNA ANALYZER - AI-POWERED VERSION 🤖\n")
    
    # Test with a sample sequence
    print("="*70)
    print("DEMO: Complete Analysis with AI Insights")
    print("="*70)
    
    test_sequence = "ATCGATCGATCGGGCCCCAAAATTTTGCGCTAGCTAGCTAGCTACGTACGTACG"
    generate_dna_report("Demo_Sequence", test_sequence, save_report=True)
