"""
FractoMap - Results Page
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent / "src"))
from common import init_session_state, get_activity_color

# Initialize
init_session_state()

# Page setup
st.set_page_config(page_title="FractoMap - Results", page_icon="📊", layout="wide")

st.markdown("""
<style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .metric-card {
        background: #f8f9fa;
        padding: 1.5rem;
        border-radius: 10px;
        text-align: center;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .metric-value {
        font-size: 2rem;
        font-weight: 700;
        color: #1E3A5F;
    }
    .metric-label {
        font-size: 0.9rem;
        color: #666;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("# 📊 Results")
st.markdown("View calculated inhibition data and statistics.")

st.markdown("---")

if st.session_state.inhibition_df is None:
    st.warning("⚠️ No results available. Please go to **Upload Data** and click **Calculate Inhibition**.")
    st.stop()

df = st.session_state.inhibition_df
control_avg = st.session_state.control_avg

# Summary statistics
st.markdown("## 📈 Summary Statistics")

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        "Control Average",
        f"{control_avg:.4f}",
        help="Average absorbance of control wells (87-96)"
    )

with col2:
    max_inh = df['% Inhibition'].max()
    st.metric(
        "Max Inhibition",
        f"{max_inh:.1f}%",
        help="Highest inhibition observed"
    )

with col3:
    active_count = len(df[df['% Inhibition'] > 50])
    st.metric(
        "Active Fractions",
        active_count,
        help="Fractions with >50% inhibition"
    )

with col4:
    best_frac = df.loc[df['% Inhibition'].idxmax(), 'Fraction']
    st.metric(
        "Best Fraction",
        f"F{best_frac}",
        help="Fraction with highest inhibition"
    )

st.markdown("---")

# Activity breakdown
st.markdown("## 🎯 Activity Distribution")

col1, col2 = st.columns([1, 2])

with col1:
    activity_counts = df['Activity'].value_counts()
    
    fig_pie = go.Figure(data=[go.Pie(
        labels=activity_counts.index,
        values=activity_counts.values,
        hole=0.4,
        marker_colors=['#1D9E75', '#EF9F27', '#F5A623', '#E24B4A']
    )])
    
    fig_pie.update_layout(
        title="Activity Distribution",
        height=300,
        showlegend=True
    )
    
    st.plotly_chart(fig_pie, use_container_width=True)

with col2:
    # Activity summary table
    activity_summary = pd.DataFrame({
        'Activity': ['Strong (>75%)', 'Moderate (50-75%)', 'Weak (25-50%)', 'Inactive (<25%)'],
        'Count': [
            len(df[df['% Inhibition'] > 75]),
            len(df[(df['% Inhibition'] > 50) & (df['% Inhibition'] <= 75)]),
            len(df[(df['% Inhibition'] > 25) & (df['% Inhibition'] <= 50)]),
            len(df[df['% Inhibition'] <= 25])
        ],
        'Color': ['🟢', '🟡', '🟠', '🔴']
    })
    
    st.dataframe(activity_summary, use_container_width=True, hide_index=True)

st.markdown("---")

# Bar chart
st.markdown("## 📊 Inhibition by Fraction")

# Create bar chart
colors = [get_activity_color(inh) for inh in df['% Inhibition']]

fig_bar = go.Figure(data=[
    go.Bar(
        x=df['Fraction'],
        y=df['% Inhibition'],
        marker_color=colors,
        hovertemplate='Fraction: F%{x}<br>Inhibition: %{y:.1f}%<extra></extra>'
    )
])

fig_bar.update_layout(
    title="% Inhibition by Fraction",
    xaxis_title="Fraction",
    yaxis_title="% Inhibition",
    yaxis_range=[0, 100],
    height=400,
    showlegend=False
)

# Add threshold line
fig_bar.add_hline(y=50, line_dash="dash", line_color="gray", 
                   annotation_text="50% threshold", annotation_position="right")

st.plotly_chart(fig_bar, use_container_width=True)

st.markdown("---")

# Detailed results table
st.markdown("## 📋 Detailed Results")

# Add color styling
def highlight_activity(row):
    if row['% Inhibition'] > 75:
        return ['background-color: rgba(29, 158, 117, 0.2)'] * len(row)
    elif row['% Inhibition'] > 50:
        return ['background-color: rgba(239, 159, 39, 0.2)'] * len(row)
    else:
        return [''] * len(row)

# Filter options
col1, col2 = st.columns([1, 3])

with col1:
    filter_activity = st.selectbox(
        "Filter by Activity",
        options=['All', 'Strong', 'Moderate', 'Weak', 'Inactive']
    )

# Apply filter
if filter_activity != 'All':
    filtered_df = df[df['Activity'] == filter_activity]
else:
    filtered_df = df

# Display table
display_df = filtered_df[['Fraction', 'Well', 'RT (min)', 'Absorbance', '% Inhibition', 'Activity']]

st.dataframe(
    display_df.style.apply(highlight_activity, axis=1).format({
        'RT (min)': '{:.2f}',
        'Absorbance': '{:.4f}',
        '% Inhibition': '{:.1f}'
    }),
    use_container_width=True,
    height=400
)

# Download buttons
st.markdown("---")

col1, col2 = st.columns(2)

with col1:
    csv = df.to_csv(index=False)
    st.download_button(
        "📥 Download Results (CSV)",
        csv,
        "inhibition_results.csv",
        "text/csv",
        use_container_width=True
    )

with col2:
    # Active fractions only
    active_df = df[df['% Inhibition'] > 50]
    csv_active = active_df.to_csv(index=False)
    st.download_button(
        "📥 Download Active Fractions Only",
        csv_active,
        "active_fractions.csv",
        "text/csv",
        use_container_width=True
    )

# Top fractions summary
st.markdown("---")
st.markdown("## 🏆 Top 10 Active Fractions")

top_10 = df.nlargest(10, '% Inhibition')[['Fraction', 'Well', 'RT (min)', '% Inhibition', 'Activity']]

st.dataframe(
    top_10.style.format({
        'RT (min)': '{:.2f}',
        '% Inhibition': '{:.1f}'
    }),
    use_container_width=True,
    hide_index=True
)
