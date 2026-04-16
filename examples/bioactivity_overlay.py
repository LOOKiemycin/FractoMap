#!/usr/bin/env python3
"""
Bioactivity-Chromatogram Overlay Tool
=====================================

A Python tool for bioactivity-guided fractionation analysis.
Overlays LC-MS/MS chromatograms with antioxidant activity data from 
microfractionation experiments.

Author: Thapanee Pruksatrakul (Visiting Scholar) Functional Metabolomics Lab, UC Riverside
Date: 15 April 2026
License: MIT

Requirements:
    pip install pymzml matplotlib pandas numpy openpyxl

Usage:
    python bioactivity_overlay.py --mzml data.mzML --inhibition inhibition.csv --output plot.png
    
    Or use as a module:
    >>> from bioactivity_overlay import BioactivityOverlay
    >>> overlay = BioactivityOverlay()
    >>> overlay.load_mzml("data.mzML")
    >>> overlay.load_inhibition_data("inhibition.csv")
    >>> overlay.plot_overlay("output.png")
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Union


class BioactivityOverlay:
    """
    A class for overlaying LC-MS/MS chromatograms with bioactivity data.
    
    Designed for microfractionation experiments where fractions are collected
    at fixed time intervals and tested for antioxidant activity (ABTS/DPPH).
    
    Attributes:
        rt (np.ndarray): Retention times from chromatogram (minutes)
        intensity (np.ndarray): Intensity values (TIC, BPC, or EIC)
        fraction_rt (np.ndarray): Retention time midpoints for each fraction
        inhibition (np.ndarray): % Inhibition values for each fraction
    """
    
    def __init__(
        self,
        collection_start: float = 1.0,
        collection_interval: float = 7/60,  # 7 seconds in minutes
        num_fractions: int = 86,
        fraction_offset: int = -1  # Offset to correct for dead volume delay
    ):
        """
        Initialize the BioactivityOverlay object.
        
        Args:
            collection_start: Start time for fraction collection (minutes)
            collection_interval: Time interval between fractions (minutes)
            num_fractions: Number of fractions collected (excluding controls)
            fraction_offset: Fraction offset to correct for dead volume delay
                            between MS detector and fraction collector.
                            Negative values shift bars earlier (most common).
                            Default: -1 (typical ~7 sec dead volume)
        """
        self.collection_start = collection_start
        self.collection_interval = collection_interval
        self.num_fractions = num_fractions
        self.fraction_offset = fraction_offset
        
        # Data containers
        self.rt: Optional[np.ndarray] = None
        self.intensity: Optional[np.ndarray] = None
        self.fraction_rt: Optional[np.ndarray] = None
        self.inhibition: Optional[np.ndarray] = None
        self.chromatogram_type: str = "TIC"
        
        # Calculate fraction retention times
        self._calculate_fraction_rt()
    
    def _calculate_fraction_rt(self) -> None:
        """Calculate the midpoint retention time for each fraction with offset correction."""
        self.fraction_rt = np.array([
            self.collection_start + (i + 0.5 + self.fraction_offset) * self.collection_interval
            for i in range(self.num_fractions)
        ])
    
    def load_mzml(
        self,
        filepath: str,
        chromatogram_type: str = "TIC",
        target_mz: Optional[float] = None,
        mz_tolerance: float = 0.01
    ) -> 'BioactivityOverlay':
        """
        Load chromatogram data from an mzML file.
        
        Args:
            filepath: Path to the mzML file
            chromatogram_type: Type of chromatogram to extract
                - "TIC": Total Ion Chromatogram (default)
                - "BPC": Base Peak Chromatogram
                - "EIC": Extracted Ion Chromatogram (requires target_mz)
            target_mz: Target m/z for EIC extraction
            mz_tolerance: m/z tolerance for EIC (Da)
        
        Returns:
            self: For method chaining
        
        Raises:
            ImportError: If pymzml is not installed
            FileNotFoundError: If mzML file doesn't exist
        """
        try:
            import pymzml
        except ImportError:
            raise ImportError("pymzml is required. Install with: pip install pymzml")
        
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"mzML file not found: {filepath}")
        
        print(f"Loading {chromatogram_type} from: {filepath.name}")
        
        run = pymzml.run.Reader(str(filepath))
        
        rt_list = []
        intensity_list = []
        
        for spectrum in run:
            # Only process MS1 spectra for chromatogram
            if spectrum.ms_level == 1:
                rt = spectrum.scan_time_in_minutes()
                
                if chromatogram_type == "TIC":
                    intensity = spectrum.TIC
                elif chromatogram_type == "BPC":
                    peaks = spectrum.peaks("centroided")
                    intensity = max(peaks[:, 1]) if len(peaks) > 0 else 0
                elif chromatogram_type == "EIC":
                    if target_mz is None:
                        raise ValueError("target_mz required for EIC extraction")
                    peaks = spectrum.peaks("centroided")
                    # Find peaks within tolerance
                    mask = np.abs(peaks[:, 0] - target_mz) <= mz_tolerance
                    intensity = np.sum(peaks[mask, 1]) if np.any(mask) else 0
                else:
                    raise ValueError(f"Unknown chromatogram_type: {chromatogram_type}")
                
                rt_list.append(rt)
                intensity_list.append(intensity)
        
        self.rt = np.array(rt_list)
        self.intensity = np.array(intensity_list)
        self.chromatogram_type = chromatogram_type
        
        print(f"  Loaded {len(self.rt)} data points")
        print(f"  RT range: {self.rt.min():.2f} - {self.rt.max():.2f} min")
        
        return self
    
    def load_inhibition_data(
        self,
        data: Union[str, np.ndarray, List[float], pd.DataFrame],
        fraction_column: Optional[str] = None,
        inhibition_column: Optional[str] = None
    ) -> 'BioactivityOverlay':
        """
        Load % inhibition data from various sources.
        
        Args:
            data: Inhibition data as:
                - str: Path to CSV or Excel file
                - np.ndarray or List: Direct array of inhibition values
                - pd.DataFrame: DataFrame with fraction and inhibition columns
            fraction_column: Column name for fraction numbers (for DataFrame/file)
            inhibition_column: Column name for inhibition values (for DataFrame/file)
        
        Returns:
            self: For method chaining
        """
        if isinstance(data, (np.ndarray, list)):
            self.inhibition = np.array(data)[:self.num_fractions]
        
        elif isinstance(data, str):
            filepath = Path(data)
            if filepath.suffix == '.csv':
                df = pd.read_csv(filepath)
            elif filepath.suffix in ['.xlsx', '.xls']:
                df = pd.read_excel(filepath)
            else:
                raise ValueError(f"Unsupported file format: {filepath.suffix}")
            
            # Auto-detect columns if not specified
            if inhibition_column is None:
                for col in ['% Inhibition', 'Inhibition', 'inhibition', '%Inhibition']:
                    if col in df.columns:
                        inhibition_column = col
                        break
            
            if inhibition_column and inhibition_column in df.columns:
                self.inhibition = df[inhibition_column].values[:self.num_fractions]
            else:
                # Assume single column of values
                self.inhibition = df.iloc[:, 0].values[:self.num_fractions]
        
        elif isinstance(data, pd.DataFrame):
            if inhibition_column and inhibition_column in data.columns:
                self.inhibition = data[inhibition_column].values[:self.num_fractions]
            else:
                self.inhibition = data.iloc[:, 0].values[:self.num_fractions]
        
        print(f"Loaded inhibition data for {len(self.inhibition)} fractions")
        print(f"  Range: {self.inhibition.min():.1f}% - {self.inhibition.max():.1f}%")
        
        return self
    
    def load_inhibition_from_plate(
        self,
        plate_data: np.ndarray,
        control_avg: float
    ) -> 'BioactivityOverlay':
        """
        Calculate % inhibition from raw absorbance plate data.
        
        Uses the simplified method where:
        - Wells 1-86 contain fractions (measured after incubation)
        - Control average is provided separately
        
        Args:
            plate_data: 8x12 numpy array of absorbance values (serpentine order)
            control_avg: Average absorbance of control wells (T₀)
        
        Returns:
            self: For method chaining
        """
        # Convert plate to serpentine order
        inhibition = []
        for row in range(8):
            if row % 2 == 0:  # Left to right
                for col in range(12):
                    frac_num = row * 12 + col + 1
                    if frac_num <= self.num_fractions:
                        abs_val = plate_data[row, col]
                        inh = (1 - abs_val / control_avg) * 100
                        inhibition.append(inh)
            else:  # Right to left
                for col in range(11, -1, -1):
                    frac_num = row * 12 + (11 - col) + 1
                    if frac_num <= self.num_fractions:
                        abs_val = plate_data[row, col]
                        inh = (1 - abs_val / control_avg) * 100
                        inhibition.append(inh)
        
        self.inhibition = np.array(inhibition)
        print(f"Calculated inhibition for {len(self.inhibition)} fractions")
        
        return self
    
    def plot_overlay(
        self,
        output_path: Optional[str] = None,
        title: Optional[str] = None,
        figsize: Tuple[int, int] = (14, 6),
        rt_range: Optional[Tuple[float, float]] = None,
        show_fraction_lines: bool = False,
        highlight_active: bool = True,
        activity_threshold: float = 50.0,
        style: str = "default"
    ) -> plt.Figure:
        """
        Create an overlay plot of chromatogram and bioactivity data.
        
        Args:
            output_path: Path to save the figure (optional)
            title: Plot title
            figsize: Figure size in inches (width, height)
            rt_range: Retention time range to display (min, max)
            show_fraction_lines: Show vertical lines at fraction boundaries
            highlight_active: Highlight active fractions (>threshold)
            activity_threshold: % Inhibition threshold for "active" fractions
            style: Plot style ("default", "publication", "poster")
        
        Returns:
            matplotlib.figure.Figure: The generated figure
        """
        if self.rt is None or self.intensity is None:
            raise ValueError("No chromatogram data loaded. Call load_mzml() first.")
        if self.inhibition is None:
            raise ValueError("No inhibition data loaded. Call load_inhibition_data() first.")
        
        # Apply style settings
        if style == "publication":
            plt.rcParams.update({
                'font.size': 10,
                'axes.labelsize': 12,
                'axes.titlesize': 14,
                'legend.fontsize': 10
            })
        elif style == "poster":
            plt.rcParams.update({
                'font.size': 14,
                'axes.labelsize': 16,
                'axes.titlesize': 18,
                'legend.fontsize': 14
            })
        
        # Create figure with two y-axes
        fig, ax1 = plt.subplots(figsize=figsize)
        ax2 = ax1.twinx()
        
        # Set RT range
        if rt_range:
            rt_mask = (self.rt >= rt_range[0]) & (self.rt <= rt_range[1])
            rt_plot = self.rt[rt_mask]
            intensity_plot = self.intensity[rt_mask]
            
            frac_mask = (self.fraction_rt >= rt_range[0]) & (self.fraction_rt <= rt_range[1])
            frac_rt_plot = self.fraction_rt[frac_mask]
            inhibition_plot = self.inhibition[frac_mask]
        else:
            rt_plot = self.rt
            intensity_plot = self.intensity
            frac_rt_plot = self.fraction_rt
            inhibition_plot = self.inhibition
        
        # Plot chromatogram (primary y-axis)
        ax1.plot(rt_plot, intensity_plot, 'b-', linewidth=0.8, alpha=0.8, label=self.chromatogram_type)
        ax1.fill_between(rt_plot, intensity_plot, alpha=0.2, color='blue')
        ax1.set_xlabel('Retention Time (min)', fontsize=12)
        ax1.set_ylabel(f'{self.chromatogram_type} Intensity', color='blue', fontsize=12)
        ax1.tick_params(axis='y', labelcolor='blue')
        ax1.set_xlim(rt_plot.min(), rt_plot.max())
        
        # Plot inhibition data as bar chart (secondary y-axis)
        bar_width = self.collection_interval * 0.8
        
        # Color bars by activity level
        colors = []
        for inh in inhibition_plot:
            if inh > 75:
                colors.append('#28a745')  # Strong - green
            elif inh > 50:
                colors.append('#ffc107')  # Moderate - yellow
            elif inh > 25:
                colors.append('#fd7e14')  # Weak - orange
            else:
                colors.append('#dc3545')  # Inactive - red
        
        bars = ax2.bar(frac_rt_plot, inhibition_plot, width=bar_width, 
                       alpha=0.6, color=colors, edgecolor='none')
        
        ax2.set_ylabel('% Inhibition', color='#28a745', fontsize=12)
        ax2.tick_params(axis='y', labelcolor='#28a745')
        ax2.set_ylim(0, max(100, inhibition_plot.max() * 1.1))
        
        # Add horizontal threshold line
        if highlight_active:
            ax2.axhline(y=activity_threshold, color='red', linestyle='--', 
                       linewidth=1, alpha=0.7, label=f'{activity_threshold}% threshold')
        
        # Show fraction boundary lines
        if show_fraction_lines:
            for i, rt in enumerate(frac_rt_plot):
                ax1.axvline(x=rt - self.collection_interval/2, 
                           color='gray', linestyle=':', linewidth=0.3, alpha=0.5)
        
        # Highlight and annotate active fractions
        if highlight_active:
            active_fracs = np.where(inhibition_plot > activity_threshold)[0]
            for idx in active_fracs:
                frac_num = idx + 1  # 1-indexed
                ax2.annotate(f'F{frac_num}', 
                           xy=(frac_rt_plot[idx], inhibition_plot[idx]),
                           xytext=(0, 5), textcoords='offset points',
                           ha='center', va='bottom', fontsize=8,
                           color='darkgreen', fontweight='bold')
        
        # Title
        if title:
            plt.title(title, fontsize=14, fontweight='bold')
        else:
            plt.title(f'Bioactivity-Chromatogram Overlay ({self.chromatogram_type})', 
                     fontsize=14, fontweight='bold')
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='blue', alpha=0.3, label=self.chromatogram_type),
            Patch(facecolor='#28a745', alpha=0.6, label='Strong (>75%)'),
            Patch(facecolor='#ffc107', alpha=0.6, label='Moderate (50-75%)'),
            Patch(facecolor='#fd7e14', alpha=0.6, label='Weak (25-50%)'),
            Patch(facecolor='#dc3545', alpha=0.6, label='Inactive (<25%)')
        ]
        ax1.legend(handles=legend_elements, loc='upper right', framealpha=0.9)
        
        plt.tight_layout()
        
        # Save figure
        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {output_path}")
        
        return fig
    
    def plot_eic_overlay(
        self,
        mzml_path: str,
        target_mz_list: List[Tuple[float, str]],
        output_path: Optional[str] = None,
        mz_tolerance: float = 0.01,
        rt_range: Optional[Tuple[float, float]] = None,
        figsize: Tuple[int, int] = (14, 8)
    ) -> plt.Figure:
        """
        Create overlay plot with multiple EICs for compound identification.
        
        Args:
            mzml_path: Path to mzML file
            target_mz_list: List of (m/z, compound_name) tuples
            output_path: Path to save figure
            mz_tolerance: m/z tolerance for EIC extraction
            rt_range: Retention time range to display
            figsize: Figure size
        
        Returns:
            matplotlib.figure.Figure
        """
        try:
            import pymzml
        except ImportError:
            raise ImportError("pymzml required")
        
        run = pymzml.run.Reader(mzml_path)
        
        # Extract all spectra first
        spectra_data = []
        for spectrum in run:
            if spectrum.ms_level == 1:
                spectra_data.append({
                    'rt': spectrum.scan_time_in_minutes(),
                    'peaks': spectrum.peaks("centroided")
                })
        
        # Create figure
        fig, axes = plt.subplots(len(target_mz_list) + 1, 1, figsize=figsize, 
                                 sharex=True, gridspec_kw={'hspace': 0.1})
        
        # Plot inhibition on top
        ax_inh = axes[0]
        bar_width = self.collection_interval * 0.8
        
        if rt_range:
            frac_mask = (self.fraction_rt >= rt_range[0]) & (self.fraction_rt <= rt_range[1])
            frac_rt_plot = self.fraction_rt[frac_mask]
            inhibition_plot = self.inhibition[frac_mask]
        else:
            frac_rt_plot = self.fraction_rt
            inhibition_plot = self.inhibition
        
        colors = ['#28a745' if i > 50 else '#ffc107' if i > 25 else '#dc3545' 
                  for i in inhibition_plot]
        ax_inh.bar(frac_rt_plot, inhibition_plot, width=bar_width, 
                   color=colors, alpha=0.7, edgecolor='none')
        ax_inh.set_ylabel('% Inhibition')
        ax_inh.set_title('Bioactivity Profile', fontsize=10, loc='left')
        ax_inh.axhline(y=50, color='red', linestyle='--', linewidth=0.8, alpha=0.7)
        
        # Plot EICs
        color_cycle = plt.cm.tab10(np.linspace(0, 1, len(target_mz_list)))
        
        for idx, (target_mz, compound_name) in enumerate(target_mz_list):
            ax = axes[idx + 1]
            
            rt_list = []
            intensity_list = []
            
            for spec in spectra_data:
                rt_list.append(spec['rt'])
                peaks = spec['peaks']
                mask = np.abs(peaks[:, 0] - target_mz) <= mz_tolerance
                intensity = np.sum(peaks[mask, 1]) if np.any(mask) else 0
                intensity_list.append(intensity)
            
            rt_arr = np.array(rt_list)
            int_arr = np.array(intensity_list)
            
            if rt_range:
                mask = (rt_arr >= rt_range[0]) & (rt_arr <= rt_range[1])
                rt_arr = rt_arr[mask]
                int_arr = int_arr[mask]
            
            ax.plot(rt_arr, int_arr, color=color_cycle[idx], linewidth=0.8)
            ax.fill_between(rt_arr, int_arr, alpha=0.3, color=color_cycle[idx])
            ax.set_ylabel('Intensity')
            ax.set_title(f'{compound_name} (m/z {target_mz:.4f})', fontsize=10, loc='left')
        
        axes[-1].set_xlabel('Retention Time (min)')
        
        plt.tight_layout()
        
        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {output_path}")
        
        return fig
    
    def find_active_compounds(
        self,
        mzml_path: str,
        activity_threshold: float = 50.0,
        rt_tolerance: float = 0.1
    ) -> pd.DataFrame:
        """
        Identify compounds in active fractions from MS data.
        
        Args:
            mzml_path: Path to mzML file
            activity_threshold: % Inhibition threshold
            rt_tolerance: RT tolerance for matching (minutes)
        
        Returns:
            pd.DataFrame: Table of compounds in active fractions
        """
        try:
            import pymzml
        except ImportError:
            raise ImportError("pymzml required")
        
        # Find active fractions
        active_fracs = np.where(self.inhibition > activity_threshold)[0]
        
        results = []
        run = pymzml.run.Reader(mzml_path)
        
        for spectrum in run:
            if spectrum.ms_level == 1:
                rt = spectrum.scan_time_in_minutes()
                
                # Check if RT falls within any active fraction
                for frac_idx in active_fracs:
                    frac_start = self.collection_start + frac_idx * self.collection_interval
                    frac_end = frac_start + self.collection_interval
                    
                    if frac_start - rt_tolerance <= rt <= frac_end + rt_tolerance:
                        peaks = spectrum.peaks("centroided")
                        
                        # Get top 10 peaks
                        if len(peaks) > 0:
                            top_indices = np.argsort(peaks[:, 1])[-10:]
                            for i in top_indices:
                                results.append({
                                    'Fraction': frac_idx + 1,
                                    'RT (min)': round(rt, 3),
                                    'm/z': round(peaks[i, 0], 4),
                                    'Intensity': peaks[i, 1],
                                    '% Inhibition': round(self.inhibition[frac_idx], 1)
                                })
                        break
        
        df = pd.DataFrame(results)
        if len(df) > 0:
            df = df.sort_values(['Fraction', 'Intensity'], ascending=[True, False])
            df = df.drop_duplicates(subset=['Fraction', 'm/z'], keep='first')
        
        return df
    
    def export_results(
        self,
        output_path: str,
        include_rt: bool = True
    ) -> None:
        """
        Export fraction-RT-inhibition mapping to CSV/Excel.
        
        Args:
            output_path: Output file path (.csv or .xlsx)
            include_rt: Include retention time information
        """
        data = {
            'Fraction': np.arange(1, self.num_fractions + 1),
            '% Inhibition': self.inhibition
        }
        
        if include_rt:
            data['RT_start (min)'] = self.collection_start + np.arange(self.num_fractions) * self.collection_interval
            data['RT_end (min)'] = data['RT_start (min)'] + self.collection_interval
            data['RT_mid (min)'] = self.fraction_rt
        
        df = pd.DataFrame(data)
        
        # Add activity classification
        df['Activity'] = pd.cut(
            df['% Inhibition'],
            bins=[-np.inf, 25, 50, 75, np.inf],
            labels=['Inactive', 'Weak', 'Moderate', 'Strong']
        )
        
        if output_path.endswith('.xlsx'):
            df.to_excel(output_path, index=False)
        else:
            df.to_csv(output_path, index=False)
        
        print(f"Results exported to: {output_path}")


def main():
    """Command-line interface for the bioactivity overlay tool."""
    parser = argparse.ArgumentParser(
        description='Overlay LC-MS/MS chromatogram with bioactivity data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --mzml data.mzML --inhibition data.csv --output plot.png
  %(prog)s --mzml data.mzML --inhibition data.csv --rt-range 1 11 --threshold 50
        """
    )
    
    parser.add_argument('--mzml', required=True, help='Path to mzML file')
    parser.add_argument('--inhibition', required=True, help='Path to inhibition data (CSV/Excel)')
    parser.add_argument('--output', default='overlay.png', help='Output figure path')
    parser.add_argument('--chromatogram', choices=['TIC', 'BPC'], default='TIC',
                        help='Chromatogram type (default: TIC)')
    parser.add_argument('--rt-range', nargs=2, type=float, metavar=('MIN', 'MAX'),
                        help='Retention time range to display')
    parser.add_argument('--threshold', type=float, default=50.0,
                        help='Activity threshold for highlighting (default: 50)')
    parser.add_argument('--collection-start', type=float, default=1.0,
                        help='Fraction collection start time in minutes (default: 1.0)')
    parser.add_argument('--collection-interval', type=float, default=7/60,
                        help='Collection interval in minutes (default: 7/60 = 7 seconds)')
    parser.add_argument('--num-fractions', type=int, default=86,
                        help='Number of fractions (default: 86)')
    parser.add_argument('--fraction-offset', type=int, default=-1,
                        help='Fraction offset to correct for dead volume delay (default: -1)')
    parser.add_argument('--title', help='Plot title')
    parser.add_argument('--style', choices=['default', 'publication', 'poster'],
                        default='default', help='Plot style')
    
    args = parser.parse_args()
    
    # Create overlay object
    overlay = BioactivityOverlay(
        collection_start=args.collection_start,
        collection_interval=args.collection_interval,
        num_fractions=args.num_fractions,
        fraction_offset=args.fraction_offset
    )
    
    # Load data
    overlay.load_mzml(args.mzml, chromatogram_type=args.chromatogram)
    overlay.load_inhibition_data(args.inhibition)
    
    # Create plot
    rt_range = tuple(args.rt_range) if args.rt_range else None
    overlay.plot_overlay(
        output_path=args.output,
        title=args.title,
        rt_range=rt_range,
        activity_threshold=args.threshold,
        style=args.style
    )
    
    print("Done!")


if __name__ == '__main__':
    main()
