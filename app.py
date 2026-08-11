"""
🧬 Ancient DNA Analyzer — Streamlit Web Dashboard
====================================================

A premium multi-page web interface for ancient DNA analysis
with interactive visualizations and AI-powered insights.

Run:  streamlit run app.py
"""

import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import streamlit as st

from src.ancient_dna.core.validator import validate_dna_sequence
from src.ancient_dna.core.stats import calculate_dna_stats
from src.ancient_dna.core.fasta_parser import read_fasta_file
from src.ancient_dna.core.transcription import transcribe_dna_to_rna, reverse_complement
from src.ancient_dna.genomics.codon_analysis import codon_usage_table, codon_bias_score
from src.ancient_dna.genomics.orf_finder import find_orfs, orf_summary
from src.ancient_dna.genomics.mutation_detector import detect_mutations, mutation_rate
from src.ancient_dna.genomics.restriction_enzymes import restriction_map
from src.ancient_dna.genomics.protein_translation import six_frame_translation, find_proteins
from src.ancient_dna.genomics.gc_analysis import gc_skew, sliding_window_gc, cumulative_gc_skew
from src.ancient_dna.visualization.plots import (
    plot_base_composition_radar,
    plot_gc_gauge,
    plot_nucleotide_heatmap,
    plot_codon_usage_bubble,
    plot_orf_map,
    plot_gc_skew_curve,
    plot_mutation_comparison,
    plot_restriction_map,
    create_full_dashboard,
)
from src.ancient_dna.ai.insights import generate_ai_insights

# =============================================================================
# PAGE CONFIG
# =============================================================================
st.set_page_config(
    page_title="Ancient DNA Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# =============================================================================
# CUSTOM CSS
# =============================================================================
st.markdown("""
<style>
    /* Dark premium theme */
    .stApp {
        background: linear-gradient(135deg, #0F172A 0%, #1E293B 50%, #0F172A 100%);
    }
    .main .block-container {
        padding-top: 2rem;
        max-width: 1200px;
    }

    /* Metric cards */
    div[data-testid="stMetric"] {
        background: linear-gradient(135deg, #1E293B, #334155);
        border: 1px solid #475569;
        border-radius: 12px;
        padding: 16px;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
    }
    div[data-testid="stMetric"] label {
        color: #94A3B8 !important;
    }
    div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
        color: #F1F5F9 !important;
        font-size: 1.8rem !important;
    }

    /* Sidebar */
    section[data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1E293B, #0F172A) !important;
        border-right: 1px solid #334155;
    }

    /* Headers */
    h1, h2, h3 { color: #F1F5F9 !important; }

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        background-color: #1E293B;
        border-radius: 8px;
        color: #94A3B8;
        border: 1px solid #334155;
        padding: 8px 20px;
    }
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #3B82F6, #8B5CF6) !important;
        color: white !important;
        border: none !important;
    }

    /* Expander */
    .streamlit-expanderHeader {
        background-color: #1E293B !important;
        border: 1px solid #334155 !important;
        border-radius: 8px !important;
        color: #F1F5F9 !important;
    }

    /* Code blocks */
    .stCodeBlock {
        border: 1px solid #334155;
        border-radius: 8px;
    }

    /* Success/info/warning */
    .stAlert {
        border-radius: 8px;
    }
</style>
""", unsafe_allow_html=True)

# =============================================================================
# SIDEBAR
# =============================================================================
with st.sidebar:
    st.markdown("## 🧬 Ancient DNA Analyzer")
    st.markdown("##### AI-Powered Paleogenomics Toolkit")
    st.divider()

    page = st.radio(
        "Navigation",
        ["🏠 Dashboard", "🔬 Genome Engineering", "🔀 Sequence Comparison", "🤖 AI Insights"],
        label_visibility="collapsed",
    )

    st.divider()
    st.markdown("### 📁 Input")

    input_method = st.radio("Input method", ["Paste sequence", "Upload FASTA"], label_visibility="collapsed")

    sequences = {}
    if input_method == "Paste sequence":
        seq_name = st.text_input("Sequence name", value="My_Sequence")
        seq_input = st.text_area("DNA Sequence", height=150, placeholder="ATCGATCG...")
        if seq_input:
            sequences[seq_name] = seq_input.upper().replace('\n', '').replace(' ', '')
    else:
        uploaded = st.file_uploader("Upload FASTA", type=['fasta', 'fa', 'fna', 'txt'])
        if uploaded:
            content = uploaded.read().decode('utf-8')
            # Parse FASTA from string
            current_name = None
            current_seq = []
            for line in content.split('\n'):
                line = line.strip()
                if line.startswith('>'):
                    if current_name:
                        sequences[current_name] = ''.join(current_seq)
                    current_name = line[1:].strip()
                    current_seq = []
                elif line and current_name:
                    current_seq.append(line.upper())
            if current_name:
                sequences[current_name] = ''.join(current_seq)

    # Demo data button
    st.divider()
    if st.button("🎯 Load Demo Data", use_container_width=True):
        sequences = {
            "Neanderthal_mtDNA": (
                "AACTCAAAGAAACCTCCTCACTATTAATCCCCATACTATTAATCATCATACTA"
                "GGCCTCGAATCAACAAATCCAATGGCACTCAACCTTCAACTAGTAATATTCC"
                "TGCCAGGACTATTCCTAACAATACTAAGCTCAGGCTGAGCCTCAAACTCAAA"
                "ATACGCCTTAATTGGAGCCCTCCGAGCAGTAGCCCAAACAATTTCATATGAAG"
                "TAAGCCTAGGTCTAATTATTCTAAGCACTATTATATTCACAGGAGGCTTCACC"
                "CTCTCAACATTTAACACCACACAAGAAACAATCTGACTAATCTTCCCAGCCTG"
                "ACCACTGGCCATAATATGATACATCTCAACCCTAGCAGAAACCAACCGAGCTC"
                "CCTTCGACCTTACAGAAGGAGAATCAGAACTTGTCTCAGGATTCAACGTCGAA"
            ),
            "Mammoth_CYTB": (
                "ATGACCAACATTCGAAAATCCCACCCACTAATAAAAATCGCTAACGATGCACT"
                "AGTCGATCTCCCCACACCCTCCAACATCTCAGCATGATGAAACTTTGGATCAC"
                "TCCTAGGCCTTTGCTTAATTACACAAATTCTAACAGGATTATTCCTAGCCATG"
                "CACTATACATCCGACATCTCCATAGCCTTCTCATCAGTAGCACACATCTGCCGA"
                "GACGTAAATTATGGTTGACTCATCCGAAACATACATGCCAACGGAGCATCCTT"
                "CTTCTTCATCTGTATCTATTTACACATTGGACGAGGCCTATATTACGGCTCTT"
                "ACCTCTATAAAGAAACCTGAAACATCGGCGTAATTCTTCTA"
            ),
        }
        st.success("Demo data loaded!")

    # Show loaded sequences
    if sequences:
        st.divider()
        st.markdown("### 📋 Loaded Sequences")
        for name, seq in sequences.items():
            st.markdown(f"**{name}** — {len(seq)} bp")

    st.divider()
    st.markdown(
        '<p style="color:#64748B; font-size:12px;">v2.0 • Built with Streamlit & Plotly</p>',
        unsafe_allow_html=True,
    )


# =============================================================================
# HELPER: Get first sequence
# =============================================================================
def get_first_seq():
    if not sequences:
        return None, None
    name = list(sequences.keys())[0]
    return name, sequences[name]


# =============================================================================
# PAGE: DASHBOARD
# =============================================================================
if page == "🏠 Dashboard":
    st.markdown("# 🏠 Analysis Dashboard")

    if not sequences:
        st.info("👈 Paste a DNA sequence or upload a FASTA file in the sidebar to get started. Or click **Load Demo Data**!")
        st.stop()

    # Analyze each sequence
    for seq_name, seq in sequences.items():
        try:
            validate_dna_sequence(seq)
        except ValueError as e:
            st.error(f"❌ {seq_name}: {e}")
            continue

        stats = calculate_dna_stats(seq)
        if not stats:
            continue

        st.markdown(f"## 📊 {seq_name}")

        # Metric cards row
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Length", f"{stats['length']} bp")
        c2.metric("GC Content", f"{stats['gc_content']}%", delta=f"{stats['gc_content'] - 50:.1f}% from 50%")
        c3.metric("AT/GC Ratio", f"{stats['at_gc_ratio']}")
        c4.metric("Pu/Py Ratio", f"{stats['purine_pyrimidine_ratio']}")
        orfs = orf_summary(seq)
        c5.metric("ORFs", f"{orfs['total_orfs']}")

        # Charts
        tab1, tab2, tab3, tab4 = st.tabs(["🧬 Composition", "🔬 GC Analysis", "🗺️ ORF Map", "✂️ Restriction"])

        with tab1:
            col1, col2 = st.columns(2)
            with col1:
                fig_radar = plot_base_composition_radar(stats, f"{seq_name} — Base Composition")
                st.plotly_chart(fig_radar, use_container_width=True)
            with col2:
                fig_gauge = plot_gc_gauge(stats['gc_content'], f"{seq_name} — GC Content")
                st.plotly_chart(fig_gauge, use_container_width=True)

            fig_heatmap = plot_nucleotide_heatmap(seq, title=f"{seq_name} — Nucleotide Distribution")
            st.plotly_chart(fig_heatmap, use_container_width=True)

        with tab2:
            skew_data = gc_skew(seq, window_size=min(50, len(seq) // 3), step_size=max(1, len(seq) // 100))
            cum_data = cumulative_gc_skew(seq)
            fig_skew = plot_gc_skew_curve(skew_data, cum_data, title=f"{seq_name} — GC Skew")
            st.plotly_chart(fig_skew, use_container_width=True)

            gc_data = sliding_window_gc(seq, window_size=min(50, len(seq) // 3))
            st.markdown("**Sliding Window GC Content (sampled)**")
            st.dataframe(
                [{"Position": d['position'], "GC %": d['gc_content']} for d in gc_data[:20]],
                use_container_width=True,
            )

        with tab3:
            orf_list = find_orfs(seq, min_length=30)
            if orf_list:
                fig_orf = plot_orf_map(orf_list, len(seq), title=f"{seq_name} — Open Reading Frames")
                st.plotly_chart(fig_orf, use_container_width=True)
                with st.expander(f"📋 ORF Details ({len(orf_list)} found)"):
                    for i, orf in enumerate(orf_list[:10]):
                        st.markdown(
                            f"**ORF {i+1}** | Frame {orf['frame']} | "
                            f"{orf['start']}–{orf['end']} bp | {orf['length']} bp "
                            f"({orf['length']//3} aa)"
                        )
            else:
                st.info("No ORFs found with minimum length 30 bp. Try adjusting parameters.")

        with tab4:
            rmap = restriction_map(seq)
            if rmap['total_cut_sites'] > 0:
                fig_rmap = plot_restriction_map(rmap, title=f"{seq_name} — Restriction Map")
                st.plotly_chart(fig_rmap, use_container_width=True)
                st.markdown(f"**{rmap['total_cut_sites']}** cut sites from **{len(rmap['enzymes_found'])}** enzymes")
                with st.expander("🔍 Enzyme Details"):
                    for enz, positions in rmap['sites'].items():
                        st.markdown(f"**{enz}**: cuts at positions {positions}")
            else:
                st.info("No restriction enzyme sites found in this sequence.")

        st.divider()


# =============================================================================
# PAGE: GENOME ENGINEERING
# =============================================================================
elif page == "🔬 Genome Engineering":
    st.markdown("# 🔬 Genome Engineering Tools")

    name, seq = get_first_seq()
    if not name:
        st.info("👈 Load a sequence first!")
        st.stop()

    st.markdown(f"**Analyzing:** {name} ({len(seq)} bp)")

    tab1, tab2, tab3, tab4 = st.tabs(["🧪 Codon Analysis", "🧬 Protein Translation", "📐 Transcription", "📊 Advanced Stats"])

    with tab1:
        st.markdown("### Codon Usage Analysis")
        frame = st.selectbox("Reading Frame", [0, 1, 2], format_func=lambda x: f"Frame {x+1}")
        table = codon_usage_table(seq, reading_frame=frame)
        bias = codon_bias_score(seq, reading_frame=frame)

        st.metric("Codon Bias Score", f"{bias}", help="20 = extreme bias, 61 = no bias")

        fig_bubble = plot_codon_usage_bubble(table, title=f"{name} — Codon Usage (Frame {frame+1})")
        st.plotly_chart(fig_bubble, use_container_width=True)

        with st.expander("📋 Full Codon Table"):
            used_codons = [c for c in table if c['count'] > 0]
            st.dataframe(used_codons, use_container_width=True)

    with tab2:
        st.markdown("### 6-Frame Protein Translation")
        translations = six_frame_translation(seq)
        for frame_label, protein in translations.items():
            with st.expander(f"Frame {frame_label} — {len(protein)} amino acids"):
                st.code(protein[:200] + ("..." if len(protein) > 200 else ""), language=None)

        st.markdown("### Potential Proteins")
        min_aa = st.slider("Minimum protein length (amino acids)", 5, 50, 10)
        proteins = find_proteins(seq, min_length=min_aa)
        if proteins:
            for i, p in enumerate(proteins[:10]):
                st.markdown(
                    f"**Protein {i+1}** | Frame {p['frame']} | "
                    f"{p['length']} aa | Start: position {p['start_aa']}"
                )
                st.code(p['sequence'][:100] + ("..." if len(p['sequence']) > 100 else ""), language=None)
        else:
            st.info(f"No proteins found with minimum length {min_aa} aa.")

    with tab3:
        st.markdown("### DNA → RNA Transcription")
        rna = transcribe_dna_to_rna(seq)
        st.code(rna[:300] + ("..." if len(rna) > 300 else ""), language=None)

        st.markdown("### Reverse Complement")
        rc = reverse_complement(seq)
        st.code(rc[:300] + ("..." if len(rc) > 300 else ""), language=None)

    with tab4:
        stats = calculate_dna_stats(seq)
        if stats:
            st.markdown("### Dinucleotide Frequencies")
            di = stats.get('dinucleotide_freq', {})
            if di:
                import plotly.graph_objects as go
                fig = go.Figure(go.Bar(
                    x=list(di.keys()), y=list(di.values()),
                    marker=dict(
                        color=list(di.values()),
                        colorscale=[[0, '#0F172A'], [0.5, '#00F0FF'], [1, '#F59E0B']],
                    ),
                ))
                fig.update_layout(
                    title="Dinucleotide Frequencies",
                    xaxis_title="Dinucleotide", yaxis_title="Frequency (%)",
                    paper_bgcolor='#0F172A', plot_bgcolor='#1E293B',
                    font=dict(color='#F1F5F9'),
                    height=400,
                )
                st.plotly_chart(fig, use_container_width=True)


# =============================================================================
# PAGE: SEQUENCE COMPARISON
# =============================================================================
elif page == "🔀 Sequence Comparison":
    st.markdown("# 🔀 Sequence Comparison")

    if len(sequences) < 2:
        st.info("Load at least **2 sequences** to compare. Try the **Demo Data** button!")
        st.stop()

    seq_names = list(sequences.keys())
    col1, col2 = st.columns(2)
    with col1:
        name1 = st.selectbox("Reference Sequence", seq_names)
    with col2:
        name2 = st.selectbox("Query Sequence", [n for n in seq_names if n != name1])

    seq1 = sequences[name1]
    seq2 = sequences[name2]
    min_len = min(len(seq1), len(seq2))

    st.markdown(f"Comparing **{min_len}** aligned positions")

    rates = mutation_rate(seq1[:min_len], seq2[:min_len])
    mutations = detect_mutations(seq1[:min_len], seq2[:min_len])

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total Mutations", rates['total_mutations'])
    c2.metric("Mutation Rate", f"{rates['mutation_rate']:.4%}")
    c3.metric("Transitions", rates['transitions'])
    c4.metric("Transversions", rates['transversions'])

    fig = plot_mutation_comparison(mutations, min_len, title=f"{name1} vs {name2}")
    st.plotly_chart(fig, use_container_width=True)

    if mutations:
        with st.expander(f"📋 Mutation Details ({len(mutations)} found)"):
            st.dataframe(mutations[:50], use_container_width=True)


# =============================================================================
# PAGE: AI INSIGHTS
# =============================================================================
elif page == "🤖 AI Insights":
    st.markdown("# 🤖 AI-Powered Insights")

    name, seq = get_first_seq()
    if not name:
        st.info("👈 Load a sequence first!")
        st.stop()

    stats = calculate_dna_stats(seq)
    if not stats:
        st.error("Failed to calculate statistics.")
        st.stop()

    st.markdown(f"**Analyzing:** {name} ({stats['length']} bp)")

    model = st.text_input("Ollama Model", value="llama3.2",
                          help="Make sure Ollama is running: `ollama serve`")

    if st.button("🚀 Generate AI Analysis", use_container_width=True, type="primary"):
        with st.spinner("Generating insights... This may take a moment."):
            insights = generate_ai_insights(name, stats, seq, model=model)
        st.markdown("---")
        st.markdown(insights)
