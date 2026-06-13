"""
FractoMap - Common functions and utilities
"""

import streamlit as st
import pandas as pd
import numpy as np
import base64
import zlib
import struct
import xml.etree.ElementTree as ET

# Constants
ROWS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

def get_serpentine_mapping():
    """Generate serpentine plate mapping for 96-well plate"""
    mapping = {}
    frac = 1
    for r in range(8):
        if r % 2 == 0:  # Left to right
            for c in range(12):
                mapping[frac] = {'row': r, 'col': c, 'well': f"{ROWS[r]}{c+1}"}
                frac += 1
        else:  # Right to left
            for c in range(11, -1, -1):
                mapping[frac] = {'row': r, 'col': c, 'well': f"{ROWS[r]}{c+1}"}
                frac += 1
    return mapping

SERPENTINE = get_serpentine_mapping()

def page_setup(title: str, icon: str = "🧪"):
    """Standard page setup with custom styling"""
    st.set_page_config(
        page_title=f"FractoMap - {title}",
        page_icon=icon,
        layout="wide"
    )
    
    # Custom CSS
    st.markdown("""
    <style>
        .main-header {
            font-size: 2rem;
            font-weight: 700;
            color: #1E3A5F;
            margin-bottom: 1rem;
        }
        .info-box {
            background: #e3f2fd;
            padding: 1rem;
            border-radius: 8px;
            border-left: 4px solid #2196f3;
            margin-bottom: 1rem;
        }
        .success-box {
            background: #e8f5e9;
            padding: 1rem;
            border-radius: 8px;
            border-left: 4px solid #4caf50;
        }
        .warning-box {
            background: #fff3e0;
            padding: 1rem;
            border-radius: 8px;
            border-left: 4px solid #ff9800;
        }
        .metric-card {
            background: #f8f9fa;
            padding: 1rem;
            border-radius: 8px;
            text-align: center;
        }
        #MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
    </style>
    """, unsafe_allow_html=True)

def decode_binary(base64_string: str, is_64bit: bool = True, is_zlib: bool = False) -> tuple:
    """Decode base64 binary data from mzML file"""
    decoded = base64.b64decode(base64_string)
    
    if is_zlib:
        try:
            decoded = zlib.decompress(decoded)
        except:
            pass
    
    if is_64bit:
        count = len(decoded) // 8
        return struct.unpack(f'<{count}d', decoded)
    else:
        count = len(decoded) // 4
        return struct.unpack(f'<{count}f', decoded)

def parse_mzml(content: str) -> dict:
    """Parse mzML file and extract TIC/BPC chromatograms"""
    try:
        root = ET.fromstring(content)
        result = {'tic': [], 'bpc': []}
        
        for elem in root.iter():
            if 'chromatogram' in elem.tag.lower():
                chrom_id = elem.get('id', '').lower()
                
                time_array = None
                intensity_array = None
                is_64bit = True
                is_zlib = False
                
                for binary_array in elem.iter():
                    if 'binaryDataArray' in binary_array.tag:
                        is_time = False
                        is_intensity = False
                        
                        for cv in binary_array.iter():
                            if 'cvParam' in cv.tag:
                                name = cv.get('name', '').lower()
                                if '64-bit' in name:
                                    is_64bit = True
                                elif '32-bit' in name:
                                    is_64bit = False
                                if 'zlib' in name:
                                    is_zlib = True
                                if 'time array' in name:
                                    is_time = True
                                if 'intensity array' in name:
                                    is_intensity = True
                        
                        for binary in binary_array.iter():
                            if 'binary' in binary.tag.lower() and binary.text:
                                data = decode_binary(binary.text.strip(), is_64bit, is_zlib)
                                if is_time:
                                    time_array = data
                                elif is_intensity:
                                    intensity_array = data
                
                if time_array and intensity_array:
                    # Convert seconds to minutes if needed
                    times = [t / 60 if t > 100 else t for t in time_array]
                    data = list(zip(times, intensity_array))
                    
                    if 'tic' in chrom_id or 'total' in chrom_id:
                        result['tic'] = data
                    elif 'bpc' in chrom_id or 'base' in chrom_id:
                        result['bpc'] = data
                    elif not result['tic']:
                        result['tic'] = data
        
        return result
    except Exception as e:
        st.error(f"Error parsing mzML: {e}")
        return {'tic': [], 'bpc': []}

def is_numeric(val):
    """Check if value is numeric"""
    if pd.isna(val):
        return False
    if isinstance(val, (int, float, np.integer, np.floating)):
        return True
    if isinstance(val, str):
        try:
            float(val)
            return True
        except:
            return False
    return False

def to_float(val):
    """Convert value to float safely"""
    if pd.isna(val):
        return 0.0
    if isinstance(val, (int, float, np.integer, np.floating)):
        return float(val)
    if isinstance(val, str):
        try:
            return float(val)
        except:
            return 0.0
    return 0.0

def parse_plate_excel(file) -> list:
    """
    Parse 96-well plate Excel file and extract 8x12 matrix.
    Supports multiple formats:
    - Simple 8x12 numeric matrix
    - Matrix with row labels (A-H) in first column
    - Matrix with header rows above data
    """
    try:
        if file.name.endswith('.csv'):
            df = pd.read_csv(file, header=None)
        else:
            df = pd.read_excel(file, header=None)
        
        matrix = []
        
        # Scan through all rows to find 8x12 numeric block
        for i in range(len(df)):
            row = df.iloc[i].values
            
            # Try to find 12 consecutive numeric values
            numeric_values = []
            
            # Method 1: Check if first column is row label (A-H)
            first_val = str(row[0]).strip().upper() if pd.notna(row[0]) else ''
            
            if first_val in ROWS:
                # Row label in first column, data in columns 1-12
                for j in range(1, min(13, len(row))):
                    if is_numeric(row[j]):
                        numeric_values.append(to_float(row[j]))
            else:
                # No row label, check first 12 columns
                for j in range(min(12, len(row))):
                    if is_numeric(row[j]):
                        numeric_values.append(to_float(row[j]))
            
            # If we found 12 numeric values, add to matrix
            if len(numeric_values) == 12:
                matrix.append(numeric_values)
                if len(matrix) == 8:
                    break
            # If we found 10-11 values (some missing), pad with zeros
            elif len(numeric_values) >= 10 and len(matrix) < 8:
                while len(numeric_values) < 12:
                    numeric_values.append(0.0)
                matrix.append(numeric_values)
                if len(matrix) == 8:
                    break
        
        if len(matrix) == 8:
            return matrix
        else:
            return None
            
    except Exception as e:
        st.error(f"Error parsing Excel: {e}")
        return None

def calculate_inhibition(plate_data: list, params: dict) -> tuple:
    """Calculate % inhibition from plate data"""
    start = params['collection_start']
    interval = params['collection_interval']
    offset = params['fraction_offset']
    
    # Calculate control average (wells 87-96)
    control_values = []
    for f in range(87, 97):
        r, c = SERPENTINE[f]['row'], SERPENTINE[f]['col']
        control_values.append(plate_data[r][c])
    control_avg = np.mean(control_values)
    
    # Calculate inhibition for each fraction
    results = []
    for f in range(1, 87):
        r, c = SERPENTINE[f]['row'], SERPENTINE[f]['col']
        absorbance = plate_data[r][c]
        inhibition = max(0, (1 - absorbance / control_avg) * 100)
        
        # Calculate retention time
        rt_center = start + (f - 1 + 0.5 + offset) * (interval / 60)
        
        # Activity classification
        if inhibition > 75:
            activity = 'Strong'
        elif inhibition > 50:
            activity = 'Moderate'
        elif inhibition > 25:
            activity = 'Weak'
        else:
            activity = 'Inactive'
        
        results.append({
            'Fraction': f,
            'Well': SERPENTINE[f]['well'],
            'Row': SERPENTINE[f]['row'],
            'Col': SERPENTINE[f]['col'],
            'RT (min)': round(rt_center, 3),
            'Absorbance': round(absorbance, 4),
            '% Inhibition': round(inhibition, 1),
            'Activity': activity
        })
    
    return pd.DataFrame(results), control_avg

def generate_demo_data() -> dict:
    """Generate demo chromatogram data"""
    tic = []
    bpc = []
    
    for rt in np.arange(0, 12, 0.02):
        # Base signal with noise
        intensity = 1e6 + np.random.random() * 1e5
        
        # Add peaks
        intensity += 8e6 * np.exp(-((rt - 5.8) / 0.15) ** 2)  # Main peak
        intensity += 3e6 * np.exp(-((rt - 4.5) / 0.2) ** 2)   # Secondary peak
        intensity += 2e6 * np.exp(-((rt - 7.2) / 0.18) ** 2)  # Third peak
        
        tic.append((rt, intensity))
        bpc.append((rt, intensity * 0.5))
    
    return {'tic': tic, 'bpc': bpc}

def get_activity_color(inhibition: float) -> str:
    """Get color based on % inhibition"""
    if inhibition > 75:
        return '#1D9E75'  # Green - Strong
    elif inhibition > 50:
        return '#EF9F27'  # Yellow - Moderate
    elif inhibition > 25:
        return '#F5A623'  # Orange - Weak
    else:
        return '#E24B4A'  # Red - Inactive

def init_session_state():
    """Initialize session state variables"""
    if 'plate_data' not in st.session_state:
        st.session_state.plate_data = None
    if 'chrom_data' not in st.session_state:
        st.session_state.chrom_data = {'tic': [], 'bpc': []}
    if 'inhibition_df' not in st.session_state:
        st.session_state.inhibition_df = None
    if 'control_avg' not in st.session_state:
        st.session_state.control_avg = None
    if 'params' not in st.session_state:
        st.session_state.params = {
            'collection_start': 1.0,
            'collection_interval': 7,
            'fraction_offset': -1,
            'assay_type': 'ABTS'
        }
