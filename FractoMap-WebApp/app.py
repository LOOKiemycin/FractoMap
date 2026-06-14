"""
FractoMap - Bioactivity-Guided Microfractionation Analysis
Single-file Streamlit App
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import base64
import zlib
import struct
import xml.etree.ElementTree as ET

# Page Config
st.set_page_config(page_title="FractoMap", page_icon="🧪", layout="wide")

# Constants
ROWS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

def get_serpentine_mapping():
    mapping = {}
    frac = 1
    for r in range(8):
        if r % 2 == 0:
            for c in range(12):
                mapping[frac] = {'row': r, 'col': c, 'well': f"{ROWS[r]}{c+1}"}
                frac += 1
        else:
            for c in range(11, -1, -1):
                mapping[frac] = {'row': r, 'col': c, 'well': f"{ROWS[r]}{c+1}"}
                frac += 1
    return mapping

SERPENTINE = get_serpentine_mapping()

# Helper Functions
def decode_binary(base64_string, is_64bit=True, is_zlib=False):
    decoded = base64.b64decode(base64_string)
    if is_zlib:
        try: decoded = zlib.decompress(decoded)
        except: pass
    if is_64bit:
        return struct.unpack(f'<{len(decoded)//8}d', decoded)
    return struct.unpack(f'<{len(decoded)//4}f', decoded)

def parse_mzml(content):
    try:
        root = ET.fromstring(content)
        result = {'tic': [], 'bpc': []}
        for elem in root.iter():
            if 'chromatogram' in elem.tag.lower():
                chrom_id = elem.get('id', '').lower()
                time_array, intensity_array = None, None
                is_64bit, is_zlib = True, False
                for ba in elem.iter():
                    if 'binaryDataArray' in ba.tag:
                        is_time, is_int = False, False
                        for cv in ba.iter():
                            if 'cvParam' in cv.tag:
                                n = cv.get('name', '').lower()
                                if '64-bit' in n: is_64bit = True
                                elif '32-bit' in n: is_64bit = False
                                if 'zlib' in n: is_zlib = True
                                if 'time array' in n: is_time = True
                                if 'intensity array' in n: is_int = True
                        for b in ba.iter():
                            if 'binary' in b.tag.lower() and b.text:
                                d = decode_binary(b.text.strip(), is_64bit, is_zlib)
                                if is_time: time_array = d
                                elif is_int: intensity_array = d
                if time_array and intensity_array:
                    times = [t/60 if t > 100 else t for t in time_array]
                    data = list(zip(times, intensity_array))
                    if 'tic' in chrom_id or 'total' in chrom_id: result['tic'] = data
                    elif 'bpc' in chrom_id or 'base' in chrom_id: result['bpc'] = data
                    elif not result['tic']: result['tic'] = data
        return result
    except Exception as e:
        st.error(f"Error parsing mzML: {e}")
        return {'tic': [], 'bpc': []}

def parse_plate_excel(file):
    try:
        df = pd.read_csv(file, header=None) if file.name.endswith('.csv') else pd.read_excel(file, header=None)
        matrix = []
        for i in range(len(df)):
            row = df.iloc[i].values
            nums = []
            start = 1 if (len(row) > 0 and isinstance(row[0], str) and row[0].strip().upper() in ROWS) else 0
            for j in range(start, min(start+12, len(row))):
                if pd.notna(row[j]):
                    try: nums.append(float(row[j]))
                    except: pass
            if len(nums) >= 10:
                while len(nums) < 12: nums.append(0.0)
                matrix.append(nums[:12])
                if len(matrix) == 8: break
        return matrix if len(matrix) == 8 else None
    except Exception as e:
        st.error(f"Error: {e}")
        return None

def calculate_inhibition(plate_data, params):
    start, interval, offset = params['start'], params['interval'], params['offset']
    ctrl = [plate_data[SERPENTINE[f]['row']][SERPENTINE[f]['col']] for f in range(87, 97)]
    ctrl_avg = np.mean(ctrl)
    results = []
    for f in range(1, 87):
        r, c = SERPENTINE[f]['row'], SERPENTINE[f]['col']
        abs_val = plate_data[r][c]
        inh = max(0, (1 - abs_val / ctrl_avg) * 100)
        rt = start + (f - 1 + 0.5 + offset) * (interval / 60)
        act = 'Strong' if inh > 75 else 'Moderate' if inh > 50 else 'Weak' if inh > 25 else 'Inactive'
        results.append({'Fraction': f, 'Well': SERPENTINE[f]['well'], 'RT (min)': round(rt, 3),
                       'Absorbance': round(abs_val, 4), '% Inhibition': round(inh, 1), 'Activity': act})
    return pd.DataFrame(results), ctrl_avg

def generate_demo():
    tic = [(rt, 1e6 + np.random.random()*1e5 + 8e6*np.exp(-((rt-5.8)/0.15)**2) + 3e6*np.exp(-((rt-4.5)/0.2)**2)) 
           for rt in np.arange(0, 12, 0.02)]
    return {'tic': tic, 'bpc': [(t, i*0.5) for t, i in tic]}

def get_color(inh):
    if inh > 75: return '#1D9E75'
    elif inh > 50: return '#EF9F27'
    elif inh > 25: return '#F5A623'
    return '#E24B4A'

# Session State
for k, v in [('plate', None), ('chrom', {'tic': [], 'bpc': []}), ('df', None), ('ctrl', None)]:
    if k not in st.session_state: st.session_state[k] = v

# Sidebar
with st.sidebar:
    st.markdown("# 🧪 FractoMap")
    st.markdown("---")
    start = st.number_input("Collection start (min)", value=1.0, step=0.1)
    interval = st.number_input("Interval (sec)", value=7, step=1)
    offset = st.number_input("Fraction offset", value=-1, step=1)
    params = {'start': start, 'interval': interval, 'offset': offset}
    st.markdown("---")
    st.markdown("[GitHub](https://github.com/LOOKiemycin/FractoMap)")

# Main
st.title("🧪 FractoMap")
tab1, tab2, tab3, tab4 = st.tabs(["📤 Upload", "🧫 Plate", "📊 Results", "📈 Overlay"])

with tab1:
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("📋 Plate Data")
        pf = st.file_uploader("Upload Excel (8×12)", type=['xlsx', 'xls', 'csv'])
        if pf:
            p = parse_plate_excel(pf)
            if p: st.session_state.plate = p; st.success("✅ Loaded")
            else: st.error("❌ Invalid format")
    with c2:
        st.subheader("📈 Chromatogram")
        cf = st.file_uploader("Upload mzML/CSV", type=['mzML', 'mzml', 'csv'])
        if cf:
            if cf.name.lower().endswith('.mzml'):
                st.session_state.chrom = parse_mzml(cf.read().decode('utf-8'))
            else:
                d = pd.read_csv(cf)
                st.session_state.chrom['tic'] = list(zip(d.iloc[:,0], d.iloc[:,1]))
            st.success(f"✅ {len(st.session_state.chrom.get('tic',[]))} points")
        if st.button("🎲 Demo Data"):
            st.session_state.chrom = generate_demo()
            st.success("✅ Demo loaded")
    
    st.markdown("---")
    if st.button("🧮 Calculate", type="primary", use_container_width=True):
        if st.session_state.plate:
            df, ctrl = calculate_inhibition(st.session_state.plate, params)
            st.session_state.df, st.session_state.ctrl = df, ctrl
            st.success("✅ Done!")
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("Control", f"{ctrl:.4f}")
            c2.metric("Max", f"{df['% Inhibition'].max():.1f}%")
            c3.metric("Active", len(df[df['% Inhibition']>50]))
            c4.metric("Best", f"F{df.loc[df['% Inhibition'].idxmax(),'Fraction']}")
        else: st.error("❌ Upload plate first")

with tab2:
    if st.session_state.plate:
        fig = go.Figure(go.Heatmap(z=np.array(st.session_state.plate), x=list(range(1,13)), y=ROWS, colorscale='RdYlGn_r'))
        fig.update_layout(title="Plate Heatmap", yaxis=dict(autorange='reversed'), height=350)
        st.plotly_chart(fig, use_container_width=True)
    else: st.warning("Upload plate first")

with tab3:
    if st.session_state.df is not None:
        df = st.session_state.df
        c1,c2,c3,c4 = st.columns(4)
        c1.metric("Control", f"{st.session_state.ctrl:.4f}")
        c2.metric("Max", f"{df['% Inhibition'].max():.1f}%")
        c3.metric("Active", len(df[df['% Inhibition']>50]))
        c4.metric("Best", f"F{df.loc[df['% Inhibition'].idxmax(),'Fraction']}")
        
        fig = go.Figure(go.Bar(x=df['Fraction'], y=df['% Inhibition'], marker_color=[get_color(i) for i in df['% Inhibition']]))
        fig.update_layout(yaxis_range=[0,100], height=300)
        fig.add_hline(y=50, line_dash="dash")
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(df, use_container_width=True, height=250)
        st.download_button("📥 CSV", df.to_csv(index=False), "results.csv")
    else: st.warning("Calculate first")

with tab4:
    c1, c2, c3 = st.columns([1,1,2])
    with c1: chrom_type = st.selectbox("Type", ["TIC", "BPC"])
    with c2: bar_width = st.slider("Bar width", 0.05, 0.3, 0.12, 0.01)
    
    df, chrom = st.session_state.df, st.session_state.chrom.get(chrom_type.lower(), [])
    if df is not None and chrom:
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        # Chromatogram - thinner line, behind bars
        fig.add_trace(go.Scatter(x=[p[0] for p in chrom], y=[p[1] for p in chrom], name=chrom_type, 
                                  line=dict(color='#2196F3', width=1), fill='tozeroy', fillcolor='rgba(33,150,243,0.15)'), secondary_y=True)
        # Bars - wider, more visible
        fig.add_trace(go.Bar(x=df['RT (min)'], y=df['% Inhibition'], marker_color=[get_color(i) for i in df['% Inhibition']], 
                             opacity=0.9, width=bar_width, name='Inhibition', marker_line_width=0), secondary_y=False)
        fig.update_layout(title=f"Overlay with {chrom_type}", xaxis_title="RT (min)", height=500, bargap=0,
                         legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99))
        fig.update_yaxes(title_text="% Inhibition", range=[0,100], secondary_y=False)
        fig.update_yaxes(title_text="Intensity", secondary_y=True)
        st.plotly_chart(fig, use_container_width=True)
        st.markdown("🔵 Chromatogram (area) | 🟢 Strong >75% | 🟡 Moderate 50-75% | 🔴 Inactive <50%")
        
        # Show active fractions
        active = df[df['% Inhibition'] > 50]
        if len(active) > 0:
            st.success(f"**Active fractions:** {', '.join([f'F{f}' for f in active['Fraction'].values])} (RT {active['RT (min)'].min():.2f}-{active['RT (min)'].max():.2f} min)")
    else: st.warning("Upload data and calculate first")

st.markdown("---")
st.markdown("**FractoMap** | [GitHub](https://github.com/LOOKiemycin/FractoMap)")
