"""
Premium DNA Visualizations
============================

Publication-quality interactive charts built with Plotly.
Uses a dark theme with neon bio-inspired accent colors,
smooth gradients, and rich hover tooltips.

Charts:
    1. Base Composition Radar Chart
    2. GC Content Gauge
    3. Nucleotide Heatmap (sliding window)
    4. Codon Usage Bubble Chart
    5. ORF Map (linear genome view)
    6. GC Skew Curve
    7. Mutation Comparison Dot Plot
    8. Restriction Enzyme Map
    9. Full Dashboard (multi-panel)
"""

import logging
from typing import Dict, List, Optional

import plotly.graph_objects as go
from plotly.subplots import make_subplots

logger = logging.getLogger(__name__)

# =============================================================================
# DESIGN SYSTEM
# =============================================================================
# Neon bio-inspired palette
COLORS = {
    'A': '#00F0FF',       # Cyan (Adenine)
    'T': '#FF6B9D',       # Pink (Thymine)
    'C': '#C084FC',       # Purple (Cytosine)
    'G': '#4ADE80',       # Green (Guanine)
    'accent1': '#F59E0B', # Amber
    'accent2': '#3B82F6', # Blue
    'accent3': '#EF4444', # Red
    'accent4': '#8B5CF6', # Violet
    'bg_dark': '#0F172A',
    'bg_card': '#1E293B',
    'bg_surface': '#334155',
    'text_primary': '#F1F5F9',
    'text_secondary': '#94A3B8',
    'grid': '#334155',
    'positive': '#4ADE80',
    'negative': '#F87171',
}

DARK_LAYOUT = dict(
    paper_bgcolor=COLORS['bg_dark'],
    plot_bgcolor=COLORS['bg_card'],
    font=dict(family='Inter, system-ui, sans-serif', color=COLORS['text_primary'], size=13),
    title_font=dict(size=20, color=COLORS['text_primary']),
    margin=dict(l=60, r=30, t=80, b=60),
)


def _apply_dark_theme(fig: go.Figure) -> None:
    """Apply the dark theme to a figure."""
    fig.update_layout(**DARK_LAYOUT)
    fig.update_xaxes(gridcolor=COLORS['grid'], zerolinecolor=COLORS['grid'])
    fig.update_yaxes(gridcolor=COLORS['grid'], zerolinecolor=COLORS['grid'])


# =============================================================================
# 1. BASE COMPOSITION RADAR CHART
# =============================================================================
def plot_base_composition_radar(stats: Dict, title: str = "Base Composition") -> go.Figure:
    """
    Create a radar/spider chart showing A, T, C, G distribution.

    Args:
        stats: Dictionary from calculate_dna_stats.
        title: Chart title.

    Returns:
        Plotly Figure.
    """
    categories = ['Adenine (A)', 'Thymine (T)', 'Cytosine (C)', 'Guanine (G)']
    values = [
        stats['percentages']['A'],
        stats['percentages']['T'],
        stats['percentages']['C'],
        stats['percentages']['G'],
    ]
    # Close the polygon
    categories_closed = categories + [categories[0]]
    values_closed = values + [values[0]]

    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(
        r=values_closed,
        theta=categories_closed,
        fill='toself',
        fillcolor='rgba(0, 240, 255, 0.15)',
        line=dict(color=COLORS['A'], width=3),
        marker=dict(size=10, color=[COLORS['A'], COLORS['T'], COLORS['C'], COLORS['G'], COLORS['A']]),
        text=[f"{v}%" for v in values_closed],
        textposition='top center',
        hovertemplate='%{theta}: %{r:.1f}%<extra></extra>',
        name='',
    ))

    fig.update_layout(
        title=dict(text=f"🧬 {title}", font=dict(size=22)),
        polar=dict(
            bgcolor=COLORS['bg_card'],
            radialaxis=dict(
                visible=True, range=[0, max(values) * 1.2],
                gridcolor=COLORS['grid'], tickfont=dict(color=COLORS['text_secondary']),
            ),
            angularaxis=dict(
                gridcolor=COLORS['grid'],
                tickfont=dict(color=COLORS['text_primary'], size=14),
            ),
        ),
        showlegend=False,
        **{k: v for k, v in DARK_LAYOUT.items() if k != 'plot_bgcolor'},
    )
    return fig


# =============================================================================
# 2. GC CONTENT GAUGE
# =============================================================================
def plot_gc_gauge(gc_content: float, title: str = "GC Content") -> go.Figure:
    """
    Create a speedometer-style gauge showing GC percentage with color zones.

    Args:
        gc_content: GC content percentage (0-100).
        title: Chart title.

    Returns:
        Plotly Figure.
    """
    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=gc_content,
        number=dict(suffix="%", font=dict(size=48, color=COLORS['text_primary'])),
        delta=dict(reference=50, valueformat=".1f", increasing_color=COLORS['positive'],
                   decreasing_color=COLORS['negative']),
        title=dict(text=f"🔬 {title}", font=dict(size=20, color=COLORS['text_primary'])),
        gauge=dict(
            axis=dict(range=[0, 100], tickwidth=2, tickcolor=COLORS['text_secondary'],
                      tickfont=dict(color=COLORS['text_secondary'])),
            bar=dict(color=COLORS['accent2'], thickness=0.3),
            bgcolor=COLORS['bg_surface'],
            borderwidth=2,
            bordercolor=COLORS['grid'],
            steps=[
                dict(range=[0, 25], color='rgba(239, 68, 68, 0.3)'),   # Low GC — AT-rich
                dict(range=[25, 40], color='rgba(245, 158, 11, 0.3)'),  # Moderate
                dict(range=[40, 60], color='rgba(74, 222, 128, 0.3)'),  # Balanced
                dict(range=[60, 75], color='rgba(245, 158, 11, 0.3)'),  # Moderate
                dict(range=[75, 100], color='rgba(239, 68, 68, 0.3)'),  # High GC
            ],
            threshold=dict(line=dict(color=COLORS['accent1'], width=4), thickness=0.8, value=gc_content),
        ),
    ))

    fig.update_layout(
        height=350,
        paper_bgcolor=COLORS['bg_dark'],
        font=dict(family='Inter, system-ui, sans-serif', color=COLORS['text_primary']),
        margin=dict(l=30, r=30, t=80, b=30),
    )
    return fig


# =============================================================================
# 3. NUCLEOTIDE HEATMAP (sliding window)
# =============================================================================
def plot_nucleotide_heatmap(
    dna_sequence: str,
    window_size: int = 50,
    title: str = "Nucleotide Distribution Heatmap",
) -> go.Figure:
    """
    Create a heatmap showing base distribution along the sequence using a sliding window.

    Args:
        dna_sequence: DNA sequence string.
        window_size: Size of the sliding window.
        title: Chart title.

    Returns:
        Plotly Figure.
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    bases = ['A', 'T', 'C', 'G']
    step = max(1, window_size // 5)

    positions = []
    matrix = {b: [] for b in bases}

    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i + window_size]
        positions.append(i + window_size // 2)
        for b in bases:
            pct = (window.count(b) / window_size) * 100
            matrix[b].append(round(pct, 1))

    z_data = [matrix[b] for b in bases]
    base_labels = ['Adenine (A)', 'Thymine (T)', 'Cytosine (C)', 'Guanine (G)']

    fig = go.Figure(go.Heatmap(
        z=z_data,
        x=positions,
        y=base_labels,
        colorscale=[
            [0, COLORS['bg_dark']],
            [0.25, '#1E3A5F'],
            [0.5, '#00F0FF'],
            [0.75, '#4ADE80'],
            [1.0, '#F59E0B'],
        ],
        colorbar=dict(
            title=dict(text='%', font=dict(color=COLORS['text_primary'])),
            tickfont=dict(color=COLORS['text_secondary']),
        ),
        hovertemplate='Position: %{x}<br>Base: %{y}<br>Content: %{z}%<extra></extra>',
    ))

    fig.update_layout(
        title=dict(text=f"🔥 {title}"),
        xaxis_title="Sequence Position (bp)",
        yaxis_title="",
    )
    _apply_dark_theme(fig)
    return fig


# =============================================================================
# 4. CODON USAGE BUBBLE CHART
# =============================================================================
def plot_codon_usage_bubble(
    codon_table: List[Dict],
    title: str = "Codon Usage",
) -> go.Figure:
    """
    Create a bubble chart with bubble sizes proportional to codon frequency.

    Args:
        codon_table: List of dicts from codon_usage_table().
        title: Chart title.

    Returns:
        Plotly Figure.
    """
    # Filter to codons that were actually used
    used = [c for c in codon_table if c['count'] > 0]
    if not used:
        fig = go.Figure()
        fig.add_annotation(text="No codons found", showarrow=False, font=dict(size=20))
        _apply_dark_theme(fig)
        return fig

    # Group by amino acid for color coding
    aa_set = sorted(set(c['amino_acid'] for c in used))
    aa_color_map = {}
    palette = ['#00F0FF', '#FF6B9D', '#C084FC', '#4ADE80', '#F59E0B', '#3B82F6',
               '#EF4444', '#8B5CF6', '#EC4899', '#14B8A6', '#F97316', '#6366F1',
               '#84CC16', '#06B6D4', '#E11D48', '#7C3AED', '#10B981', '#F43F5E',
               '#0EA5E9', '#A855F7', '#22D3EE']
    for i, aa in enumerate(aa_set):
        aa_color_map[aa] = palette[i % len(palette)]

    fig = go.Figure()

    for aa in aa_set:
        group = [c for c in used if c['amino_acid'] == aa]
        fig.add_trace(go.Scatter(
            x=[c['codon'] for c in group],
            y=[c['rscu'] for c in group],
            mode='markers+text',
            marker=dict(
                size=[max(c['count'] * 3, 8) for c in group],
                color=aa_color_map[aa],
                opacity=0.8,
                line=dict(width=1, color='rgba(255,255,255,0.3)'),
            ),
            text=[c['amino_acid'] for c in group],
            textposition='top center',
            textfont=dict(size=9, color=COLORS['text_secondary']),
            name=aa,
            hovertemplate='Codon: %{x}<br>RSCU: %{y:.2f}<br>Count: %{customdata}<extra>%{text}</extra>',
            customdata=[c['count'] for c in group],
        ))

    fig.update_layout(
        title=dict(text=f"🫧 {title}"),
        xaxis_title="Codon",
        yaxis_title="RSCU (Relative Synonymous Codon Usage)",
        xaxis=dict(tickangle=90, tickfont=dict(size=8)),
        showlegend=True,
        legend=dict(
            title=dict(text="Amino Acid", font=dict(color=COLORS['text_primary'])),
            bgcolor='rgba(30, 41, 59, 0.8)',
            font=dict(color=COLORS['text_secondary'], size=10),
        ),
        height=550,
    )
    _apply_dark_theme(fig)
    return fig


# =============================================================================
# 5. ORF MAP (linear genome view)
# =============================================================================
def plot_orf_map(orfs: List[Dict], seq_length: int, title: str = "ORF Map") -> go.Figure:
    """
    Create a linear genome map showing ORF positions across reading frames.

    Args:
        orfs: List of ORF dicts from find_orfs().
        seq_length: Total sequence length.
        title: Chart title.

    Returns:
        Plotly Figure.
    """
    fig = go.Figure()

    frame_colors = {
        1: '#00F0FF', 2: '#4ADE80', 3: '#F59E0B',
        -1: '#FF6B9D', -2: '#C084FC', -3: '#EF4444',
    }
    frame_y = {1: 3, 2: 2, 3: 1, -1: -1, -2: -2, -3: -3}

    # Draw frame lines
    for frame, y in frame_y.items():
        label = f"Frame {'+' if frame > 0 else ''}{frame}"
        fig.add_trace(go.Scatter(
            x=[0, seq_length], y=[y, y],
            mode='lines',
            line=dict(color=COLORS['grid'], width=1, dash='dot'),
            showlegend=False, hoverinfo='skip',
        ))
        fig.add_annotation(
            x=-seq_length * 0.02, y=y, text=label,
            showarrow=False, font=dict(size=11, color=frame_colors.get(frame, '#fff')),
            xanchor='right',
        )

    # Draw ORFs as thick bars
    for orf in orfs:
        y = frame_y.get(orf['frame'], 0)
        color = frame_colors.get(orf['frame'], '#fff')
        fig.add_trace(go.Scatter(
            x=[orf['start'], orf['end']],
            y=[y, y],
            mode='lines',
            line=dict(color=color, width=14),
            opacity=0.85,
            showlegend=False,
            hovertemplate=(
                f"Frame {orf['frame']}<br>"
                f"Start: {orf['start']} bp<br>"
                f"End: {orf['end']} bp<br>"
                f"Length: {orf['length']} bp ({orf['length']//3} aa)"
                "<extra></extra>"
            ),
        ))
        # ORF label
        mid = (orf['start'] + orf['end']) / 2
        fig.add_annotation(
            x=mid, y=y, text=f"{orf['length']}bp",
            showarrow=False, font=dict(size=9, color='white'),
        )

    # Separator line between forward and reverse
    fig.add_hline(y=0, line=dict(color=COLORS['accent1'], width=2, dash='dash'))
    fig.add_annotation(x=seq_length / 2, y=0.4, text="Forward ↑ | Reverse ↓",
                       showarrow=False, font=dict(size=11, color=COLORS['accent1']))

    fig.update_layout(
        title=dict(text=f"🗺️ {title}"),
        xaxis_title="Sequence Position (bp)",
        yaxis=dict(showticklabels=False, range=[-4, 4.5]),
        xaxis=dict(range=[-seq_length * 0.08, seq_length * 1.02]),
        height=400,
        showlegend=False,
    )
    _apply_dark_theme(fig)
    return fig


# =============================================================================
# 6. GC SKEW CURVE
# =============================================================================
def plot_gc_skew_curve(
    skew_data: List[Dict],
    cumulative_data: List[Dict] = None,
    title: str = "GC Skew Analysis",
) -> go.Figure:
    """
    Plot GC skew and optional cumulative GC skew curves.

    Args:
        skew_data: List of dicts from gc_skew().
        cumulative_data: Optional list from cumulative_gc_skew().
        title: Chart title.

    Returns:
        Plotly Figure.
    """
    has_cumulative = cumulative_data and len(cumulative_data) > 0
    rows = 2 if has_cumulative else 1

    fig = make_subplots(
        rows=rows, cols=1,
        subplot_titles=["GC Skew (G-C)/(G+C)"] + (["Cumulative GC Skew"] if has_cumulative else []),
        vertical_spacing=0.15,
    )

    positions = [d['position'] for d in skew_data]
    values = [d['skew_value'] for d in skew_data]

    # Color positive and negative differently
    colors = [COLORS['positive'] if v >= 0 else COLORS['negative'] for v in values]

    fig.add_trace(go.Bar(
        x=positions, y=values,
        marker=dict(color=colors, opacity=0.8),
        hovertemplate='Position: %{x}<br>GC Skew: %{y:.4f}<extra></extra>',
        name='GC Skew',
    ), row=1, col=1)

    if has_cumulative:
        cum_pos = [d['position'] for d in cumulative_data]
        cum_vals = [d['cumulative_skew'] for d in cumulative_data]

        fig.add_trace(go.Scatter(
            x=cum_pos, y=cum_vals,
            mode='lines',
            line=dict(color=COLORS['A'], width=2),
            fill='tozeroy',
            fillcolor='rgba(0, 240, 255, 0.1)',
            hovertemplate='Position: %{x}<br>Cumulative: %{y:.0f}<extra></extra>',
            name='Cumulative',
        ), row=2, col=1)

        # Mark min/max
        min_idx = min(range(len(cum_vals)), key=lambda i: cum_vals[i])
        max_idx = max(range(len(cum_vals)), key=lambda i: cum_vals[i])
        fig.add_annotation(
            x=cum_pos[min_idx], y=cum_vals[min_idx],
            text="Origin of Replication?", showarrow=True,
            arrowhead=2, arrowcolor=COLORS['positive'],
            font=dict(color=COLORS['positive'], size=11),
            row=2, col=1,
        )

    fig.update_layout(
        title=dict(text=f"📈 {title}"),
        height=500 if has_cumulative else 350,
        showlegend=False,
    )
    _apply_dark_theme(fig)
    return fig


# =============================================================================
# 7. MUTATION COMPARISON DOT PLOT
# =============================================================================
def plot_mutation_comparison(
    mutations: List[Dict],
    seq_length: int,
    title: str = "Mutation Map",
) -> go.Figure:
    """
    Create a visual map of mutations between two sequences.

    Args:
        mutations: List of mutation dicts from detect_mutations().
        seq_length: Length of the compared sequences.
        title: Chart title.

    Returns:
        Plotly Figure.
    """
    fig = go.Figure()

    if not mutations:
        fig.add_annotation(text="No mutations detected ✅", showarrow=False,
                           font=dict(size=20, color=COLORS['positive']))
        _apply_dark_theme(fig)
        return fig

    # Separate transitions and transversions
    transitions = [m for m in mutations if m['type'] == 'transition']
    transversions = [m for m in mutations if m['type'] == 'transversion']

    if transitions:
        fig.add_trace(go.Scatter(
            x=[m['position'] for m in transitions],
            y=[1] * len(transitions),
            mode='markers',
            marker=dict(
                size=10, color=COLORS['accent2'], symbol='circle',
                line=dict(width=1, color='rgba(255,255,255,0.3)'),
            ),
            text=[f"{m['ref_base']}→{m['alt_base']}" for m in transitions],
            hovertemplate='Pos: %{x}<br>%{text}<br>Type: Transition<extra></extra>',
            name=f'Transitions ({len(transitions)})',
        ))

    if transversions:
        fig.add_trace(go.Scatter(
            x=[m['position'] for m in transversions],
            y=[2] * len(transversions),
            mode='markers',
            marker=dict(
                size=12, color=COLORS['accent3'], symbol='diamond',
                line=dict(width=1, color='rgba(255,255,255,0.3)'),
            ),
            text=[f"{m['ref_base']}→{m['alt_base']}" for m in transversions],
            hovertemplate='Pos: %{x}<br>%{text}<br>Type: Transversion<extra></extra>',
            name=f'Transversions ({len(transversions)})',
        ))

    total = len(mutations)
    rate = round(total / seq_length * 100, 2) if seq_length > 0 else 0
    fig.add_annotation(
        x=0.5, y=1.12, xref='paper', yref='paper',
        text=f"Total: {total} mutations | Rate: {rate}%",
        showarrow=False, font=dict(size=13, color=COLORS['text_secondary']),
    )

    fig.update_layout(
        title=dict(text=f"🔀 {title}"),
        xaxis_title="Sequence Position (bp)",
        yaxis=dict(
            tickvals=[1, 2], ticktext=['Transitions', 'Transversions'],
            range=[0, 3],
        ),
        legend=dict(bgcolor='rgba(30,41,59,0.8)', font=dict(color=COLORS['text_secondary'])),
        height=350,
    )
    _apply_dark_theme(fig)
    return fig


# =============================================================================
# 8. RESTRICTION ENZYME MAP
# =============================================================================
def plot_restriction_map(
    restriction_data: Dict,
    title: str = "Restriction Map",
) -> go.Figure:
    """
    Create a linear map with enzyme cut sites marked.

    Args:
        restriction_data: Dict from restriction_map().
        title: Chart title.

    Returns:
        Plotly Figure.
    """
    fig = go.Figure()
    seq_len = restriction_data['sequence_length']
    sites = restriction_data.get('sites', {})

    if not sites:
        fig.add_annotation(text="No restriction sites found", showarrow=False,
                           font=dict(size=20, color=COLORS['text_secondary']))
        _apply_dark_theme(fig)
        return fig

    # Draw the DNA backbone
    fig.add_trace(go.Scatter(
        x=[0, seq_len], y=[0, 0],
        mode='lines',
        line=dict(color=COLORS['A'], width=6),
        showlegend=False, hoverinfo='skip',
    ))

    enzyme_list = list(sites.keys())
    palette = ['#00F0FF', '#FF6B9D', '#C084FC', '#4ADE80', '#F59E0B', '#3B82F6',
               '#EF4444', '#8B5CF6', '#EC4899', '#14B8A6', '#F97316', '#6366F1']

    for idx, (enzyme, positions) in enumerate(sites.items()):
        color = palette[idx % len(palette)]
        y_offset = 0.5 + (idx % 4) * 0.4

        for pos in positions:
            # Cut site marker
            fig.add_trace(go.Scatter(
                x=[pos, pos], y=[-0.2, y_offset],
                mode='lines+markers',
                line=dict(color=color, width=2, dash='dot'),
                marker=dict(size=[0, 8], color=color, symbol='triangle-down'),
                showlegend=False,
                hovertemplate=f'{enzyme}<br>Cut position: {pos} bp<extra></extra>',
            ))
            fig.add_annotation(
                x=pos, y=y_offset + 0.15, text=enzyme,
                showarrow=False, font=dict(size=9, color=color),
                textangle=-45,
            )

    fig.update_layout(
        title=dict(text=f"✂️ {title}"),
        xaxis_title="Sequence Position (bp)",
        yaxis=dict(showticklabels=False, range=[-0.5, 3]),
        height=350,
        showlegend=False,
    )
    _apply_dark_theme(fig)
    return fig


# =============================================================================
# 9. FULL DASHBOARD
# =============================================================================
def create_full_dashboard(
    sequence_name: str,
    stats: Dict,
    dna_sequence: str,
    save_path: Optional[str] = None,
) -> go.Figure:
    """
    Create a multi-panel dashboard combining key visualizations.

    Args:
        sequence_name: Name of the sequence.
        stats: Dict from calculate_dna_stats.
        dna_sequence: The DNA sequence string.
        save_path: If provided, saves the dashboard as an HTML file.

    Returns:
        Plotly Figure with subplots.
    """
    from src.ancient_dna.genomics.gc_analysis import sliding_window_gc

    fig = make_subplots(
        rows=2, cols=3,
        specs=[
            [{'type': 'polar'}, {'type': 'indicator'}, {'type': 'bar'}],
            [{'type': 'heatmap', 'colspan': 2}, None, {'type': 'bar'}],
        ],
        subplot_titles=[
            "Base Composition", "GC Content", "Base Counts",
            "Nucleotide Heatmap", "", "Dinucleotide Freq",
        ],
        vertical_spacing=0.18,
        horizontal_spacing=0.08,
    )

    bases = ['A', 'T', 'C', 'G']
    base_colors = [COLORS['A'], COLORS['T'], COLORS['C'], COLORS['G']]

    # Panel 1: Radar
    values = [stats['percentages'][b] for b in bases] + [stats['percentages']['A']]
    theta = ['A', 'T', 'C', 'G', 'A']
    fig.add_trace(go.Scatterpolar(
        r=values, theta=theta, fill='toself',
        fillcolor='rgba(0, 240, 255, 0.15)',
        line=dict(color=COLORS['A'], width=2),
        name='',
    ), row=1, col=1)

    # Panel 2: Gauge
    fig.add_trace(go.Indicator(
        mode="gauge+number",
        value=stats['gc_content'],
        number=dict(suffix="%"),
        gauge=dict(
            axis=dict(range=[0, 100]),
            bar=dict(color=COLORS['accent2']),
            steps=[
                dict(range=[0, 30], color='rgba(239, 68, 68, 0.3)'),
                dict(range=[30, 70], color='rgba(74, 222, 128, 0.3)'),
                dict(range=[70, 100], color='rgba(239, 68, 68, 0.3)'),
            ],
        ),
    ), row=1, col=2)

    # Panel 3: Bar chart
    fig.add_trace(go.Bar(
        x=bases,
        y=[stats['counts'][b] for b in bases],
        marker=dict(
            color=base_colors,
            line=dict(width=1, color='rgba(255,255,255,0.2)'),
        ),
        text=[stats['counts'][b] for b in bases],
        textposition='auto',
        textfont=dict(color='white'),
        hovertemplate='%{x}: %{y} bases<extra></extra>',
        name='',
    ), row=1, col=3)

    # Panel 4: Heatmap
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    window = max(20, len(seq) // 20)
    step = max(1, window // 5)
    positions = []
    matrix = {b: [] for b in bases}
    for i in range(0, len(seq) - window + 1, step):
        w = seq[i:i + window]
        positions.append(i + window // 2)
        for b in bases:
            matrix[b].append(round((w.count(b) / window) * 100, 1))

    fig.add_trace(go.Heatmap(
        z=[matrix[b] for b in bases],
        x=positions, y=bases,
        colorscale=[[0, COLORS['bg_dark']], [0.5, '#00F0FF'], [1, '#F59E0B']],
        showscale=False,
        hovertemplate='Pos: %{x}<br>%{y}: %{z}%<extra></extra>',
    ), row=2, col=1)

    # Panel 5: Dinucleotide frequencies
    di_freq = stats.get('dinucleotide_freq', {})
    if di_freq:
        top_di = sorted(di_freq.items(), key=lambda x: x[1], reverse=True)[:12]
        fig.add_trace(go.Bar(
            x=[d[0] for d in top_di],
            y=[d[1] for d in top_di],
            marker=dict(color=COLORS['accent4'], opacity=0.8),
            hovertemplate='%{x}: %{y}%<extra></extra>',
            name='',
        ), row=2, col=3)

    fig.update_layout(
        title=dict(
            text=f"🧬 {sequence_name} — DNA Analysis Dashboard",
            font=dict(size=24, color=COLORS['text_primary']),
        ),
        height=750,
        showlegend=False,
        **{k: v for k, v in DARK_LAYOUT.items() if k not in ('title_font', 'margin')},
        margin=dict(l=60, r=30, t=100, b=60),
    )
    fig.update_xaxes(gridcolor=COLORS['grid'], zerolinecolor=COLORS['grid'])
    fig.update_yaxes(gridcolor=COLORS['grid'], zerolinecolor=COLORS['grid'])

    if save_path:
        fig.write_html(save_path)
        logger.info("Dashboard saved to '%s'", save_path)

    return fig
