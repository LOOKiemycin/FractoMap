"""
FractoMap - Overlay Page
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent / "src"))
from common import init_session_state, get_activity_color

# Initialize
init_session_state()

# Page setup
st.set_page_config(page_title="FractoMap - Overlay", page_icon="📈", layout="wide")

st.markdown("""
<style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .legend-box {
        background: #f8f9fa;
        padding: 1rem;
        border-radius: 8px;
        margin-top: 1rem;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("# 📈 Bioactivity-Chromatogram Overlay")
st.markdown("Visualize bioactivity data overlaid with LC-MS chromatogram.")

st.markdown("---")

# Check data availability
if st.session_state.inhibition_df is None:
    st.warning("⚠️ No inhibition data. Please go to **Upload Data** and click **Calculate Inhibition**.")
    st.stop()

if not st.session_state.chrom_data.get('tic') and not st.session_state.chrom_data.get('bpc'):
    st.warning("⚠️ No chromatogram data. Please upload mzML file or generate demo data on **Upload Data** page.")
    st.stop()

df = st.session_state.inhibition_df
chrom_data = st.session_state.chrom_data

# Options
col1, col2, col3 = st.columns([1, 1, 2])

with col1:
    chrom_type = st.selectbox(
        "Chromatogram Type",
        options=['TIC', 'BPC'],
        help="TIC = Total Ion Current, BPC = Base Peak Chromatogram"
    )

with col2:
    show_labels = st.checkbox("Show fraction labels", value=False)

# Get chromatogram data
chrom = chrom_data.get(chrom_type.lower(), [])

if not chrom:
    st.warning(f"⚠️ No {chrom_type} data available. Try the other option.")
    st.stop()

# Create overlay plot
st.markdown("## Overlay Plot")

fig = make_subplots(specs=[[{"secondary_y": True}]])

# Chromatogram line
rt_vals = [p[0] for p in chrom]
int_vals = [p[1] for p in chrom]

fig.add_trace(
    go.Scatter(
        x=rt_vals,
        y=int_vals,
        name=chrom_type,
        line=dict(color='#2196F3', width=1.5),
        hovertemplate='RT: %{x:.2f} min<br>Intensity: %{y:.2e}<extra></extra>'
    ),
    secondary_y=True
)

# Inhibition bars
colors = [get_activity_color(inh) for inh in df['% Inhibition']]

fig.add_trace(
    go.Bar(
        x=df['RT (min)'],
        y=df['% Inhibition'],
        name='% Inhibition',
        marker_color=colors,
        opacity=0.7,
        width=0.08,
        hovertemplate='Fraction: F%{customdata}<br>RT: %{x:.2f} min<br>Inhibition: %{y:.1f}%<extra></extra>',
        customdata=df['Fraction']
    ),
    secondary_y=False
)

# Add labels for active fractions
if show_labels:
    active = df[df['% Inhibition'] > 50]
    for _, row in active.iterrows():
        fig.add_annotation(
            x=row['RT (min)'],
            y=row['% Inhibition'],
            text=f"F{row['Fraction']}",
            showarrow=True,
            arrowhead=2,
            arrowsize=1,
            arrowwidth=1,
            ax=0,
            ay=-30,
            font=dict(size=10)
        )

# Layout
fig.update_layout(
    title=f"Bioactivity Overlay with {chrom_type}",
    xaxis_title="Retention Time (min)",
    height=500,
    legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="right",
        x=0.99,
        bgcolor="rgba(255,255,255,0.8)"
    ),
    hovermode="x unified"
)

fig.update_yaxes(
    title_text="% Inhibition",
    range=[0, 100],
    secondary_y=False
)

fig.update_yaxes(
    title_text="Intensity",
    secondary_y=True,
    tickformat=".1e"
)

st.plotly_chart(fig, use_container_width=True)

# Legend
st.markdown("""
<div class="legend-box">
    <b>Legend:</b><br>
    🔵 <b>Blue line</b> = Chromatogram (TIC/BPC)<br>
    🟢 <b>Green bars</b> = Strong inhibition (>75%)<br>
    🟡 <b>Yellow bars</b> = Moderate inhibition (50-75%)<br>
    🟠 <b>Orange bars</b> = Weak inhibition (25-50%)<br>
    🔴 <b>Red bars</b> = Inactive (<25%)
</div>
""", unsafe_allow_html=True)

st.markdown("---")

# Active regions summary
st.markdown("## 🎯 Active Regions")

active_df = df[df['% Inhibition'] > 50].copy()

if len(active_df) > 0:
    # Find contiguous regions
    active_df = active_df.sort_values('Fraction')
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Active Fractions")
        st.dataframe(
            active_df[['Fraction', 'Well', 'RT (min)', '% Inhibition', 'Activity']].style.format({
                'RT (min)': '{:.2f}',
                '% Inhibition': '{:.1f}'
            }),
            use_container_width=True,
            hide_index=True
        )
    
    with col2:
        st.markdown("### RT Windows")
        
        # Group contiguous fractions
        groups = []
        current_group = [active_df.iloc[0]]
        
        for i in range(1, len(active_df)):
            if active_df.iloc[i]['Fraction'] - active_df.iloc[i-1]['Fraction'] == 1:
                current_group.append(active_df.iloc[i])
            else:
                groups.append(current_group)
                current_group = [active_df.iloc[i]]
        groups.append(current_group)
        
        # Display groups
        for i, group in enumerate(groups):
            start_frac = group[0]['Fraction']
            end_frac = group[-1]['Fraction']
            start_rt = group[0]['RT (min)']
            end_rt = group[-1]['RT (min)']
            max_inh = max(g['% Inhibition'] for g in group)
            
            st.markdown(f"""
            **Region {i+1}:** F{start_frac}-F{end_frac}  
            RT: {start_rt:.2f} - {end_rt:.2f} min  
            Max inhibition: {max_inh:.1f}%
            """)
else:
    st.info("No fractions with >50% inhibition found.")

st.markdown("---")

# Zoomed view
st.markdown("## 🔍 Zoom to Active Region")

col1, col2 = st.columns(2)

with col1:
    rt_min = st.number_input(
        "RT Start (min)",
        value=float(df['RT (min)'].min()),
        min_value=0.0,
        max_value=float(df['RT (min)'].max()),
        step=0.5
    )

with col2:
    rt_max = st.number_input(
        "RT End (min)",
        value=float(df['RT (min)'].max()),
        min_value=0.0,
        max_value=float(df['RT (min)'].max()) + 1,
        step=0.5
    )

# Zoomed plot
fig_zoom = make_subplots(specs=[[{"secondary_y": True}]])

# Filter data for zoom
rt_mask = [(rt >= rt_min and rt <= rt_max) for rt, _ in chrom]
chrom_zoom = [c for c, m in zip(chrom, rt_mask) if m]

df_zoom = df[(df['RT (min)'] >= rt_min) & (df['RT (min)'] <= rt_max)]

if chrom_zoom:
    fig_zoom.add_trace(
        go.Scatter(
            x=[p[0] for p in chrom_zoom],
            y=[p[1] for p in chrom_zoom],
            name=chrom_type,
            line=dict(color='#2196F3', width=2),
            fill='tozeroy',
            fillcolor='rgba(33, 150, 243, 0.1)'
        ),
        secondary_y=True
    )

if len(df_zoom) > 0:
    colors_zoom = [get_activity_color(inh) for inh in df_zoom['% Inhibition']]
    
    fig_zoom.add_trace(
        go.Bar(
            x=df_zoom['RT (min)'],
            y=df_zoom['% Inhibition'],
            name='% Inhibition',
            marker_color=colors_zoom,
            opacity=0.8,
            width=0.1,
            text=[f"F{f}" for f in df_zoom['Fraction']],
            textposition='outside',
            textfont=dict(size=10)
        ),
        secondary_y=False
    )

fig_zoom.update_layout(
    title=f"Zoomed View: {rt_min:.1f} - {rt_max:.1f} min",
    xaxis_title="Retention Time (min)",
    height=400,
    showlegend=False
)

fig_zoom.update_yaxes(title_text="% Inhibition", range=[0, 110], secondary_y=False)
fig_zoom.update_yaxes(title_text="Intensity", secondary_y=True, tickformat=".1e")

st.plotly_chart(fig_zoom, use_container_width=True)

# Export options
st.markdown("---")
st.markdown("## 📥 Export")

col1, col2 = st.columns(2)

with col1:
    # Export active regions
    if len(active_df) > 0:
        csv = active_df.to_csv(index=False)
        st.download_button(
            "📥 Download Active Regions (CSV)",
            csv,
            "active_regions.csv",
            "text/csv",
            use_container_width=True
        )

with col2:
    # Export full results with RT
    full_export = df[['Fraction', 'Well', 'RT (min)', 'Absorbance', '% Inhibition', 'Activity']]
    csv_full = full_export.to_csv(index=False)
    st.download_button(
        "📥 Download Full Results (CSV)",
        csv_full,
        "fractomap_results.csv",
        "text/csv",
        use_container_width=True
    )
