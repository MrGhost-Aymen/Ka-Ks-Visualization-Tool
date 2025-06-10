#!/usr/bin/env python3
"""
Robust Ka/Ks Visualization Tool with Flexible Column Handling
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from argparse import ArgumentParser
from scipy.cluster import hierarchy

def validate_dataframe(df):
    """Validate input DataFrame structure and content."""
    # Try common column name variations
    ka_ks_cols = ['Ka/Ks', 'Ka_Ks', 'KaKs', 'dN/dS', 'dn_ds']
    found_col = None
    
    for col in ka_ks_cols:
        if col in df.columns:
            found_col = col
            break
            
    if not found_col:
        print("Error: Could not find Ka/Ks column. Tried:")
        print(ka_ks_cols)
        print(f"Available columns: {list(df.columns)}")
        return None
    
    required_columns = ['Gene', 'Sequence1', 'Sequence2', found_col]
    missing_cols = [col for col in required_columns if col not in df.columns]
    
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        print(f"Available columns: {list(df.columns)}")
        return None
    
    try:
        df[found_col] = pd.to_numeric(df[found_col], errors='raise')
    except ValueError as e:
        print(f"Error in {found_col} values: {e}")
        return None
    
    return found_col

def preprocess_data(df, ka_ks_col, args):
    """Clean and prepare the data for visualization."""
    # Create consistent species pair names
    df['species_pair'] = df.apply(
        lambda row: f"{sorted([row['Sequence1'], row['Sequence2']])[0]} vs {sorted([row['Sequence1'], row['Sequence2']])[1]}",
        axis=1
    )
    
    # Create cleaned Ka/Ks column
    df['Ka_Ks_processed'] = df[ka_ks_col].replace([np.inf, -np.inf], np.nan)
    
    if args.log_transform:
        df['Ka_Ks_processed'] = np.log2(df['Ka_Ks_processed'])
        df.loc[df['Ka_Ks_processed'].isna(), 'Ka_Ks_processed'] = 0
    
    if 'p_value' in df.columns and 'Significant' not in df.columns:
        df['Significant'] = df['p_value'] < 0.05
    
    return df

def create_heatmap(pivoted_df, args, output_dir, suffix=""):
    """Generate heatmap visualization."""
    plt.figure(figsize=args.figsize)
    
    annot = pivoted_df.apply(lambda col: col.map(lambda x: f"{x:.2f}" if pd.notna(x) else ""))
    
    ax = sns.heatmap(
        pivoted_df,
        cmap=args.colormap,
        annot=annot if args.annotate else False,
        fmt="",
        linewidths=0.5,
        cbar_kws={'label': 'log2(Ka/Ks)' if args.log_transform else 'Ka/Ks'},
        mask=pivoted_df.isna(),
        vmin=args.vmin,
        vmax=args.vmax
    )
    
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/heatmap{suffix}.png", dpi=args.dpi)
    plt.close()

def create_clustered_heatmap(pivoted_df, args, output_dir):
    """Generate clustered heatmap."""
    try:
        clust_data = pivoted_df.fillna(0)
        row_linkage = hierarchy.linkage(clust_data, method=args.cluster_method)
        col_linkage = hierarchy.linkage(clust_data.T, method=args.cluster_method)
        
        annot = clust_data.apply(lambda col: col.map(lambda x: f"{x:.2f}" if x != 0 else ""))
        
        g = sns.clustermap(
            clust_data,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            cmap=args.colormap,
            annot=annot if args.annotate else False,
            fmt="",
            figsize=args.figsize,
            cbar_kws={'label': 'log2(Ka/Ks)' if args.log_transform else 'Ka/Ks'}
        )
        
        plt.savefig(f"{output_dir}/clustered_heatmap.png", dpi=args.dpi)
        plt.close()
    except Exception as e:
        print(f"Error in clustering: {e}")

def create_dot_plot(df, ka_ks_col, args, output_dir):
    """Create dot plot visualization."""
    plt.figure(figsize=(12, 8))
    
    y_col = 'Ka_Ks_processed' if 'Ka_Ks_processed' in df else ka_ks_col
    
    sns.stripplot(
        data=df,
        x='species_pair',
        y=y_col,
        hue='Gene',
        jitter=True,
        dodge=True,
        palette='viridis'
    )
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/dot_plot.png", dpi=args.dpi)
    plt.close()

def main():
    parser = ArgumentParser(description="Ka/Ks Visualization Tool")
    parser.add_argument("input_file", help="Input CSV file")
    parser.add_argument("--output_dir", default="results", help="Output directory")
    parser.add_argument("--cluster", action="store_true", help="Create clustered heatmap")
    parser.add_argument("--log_transform", action="store_true", help="Apply log2 transformation")
    parser.add_argument("--annotate", action="store_true", help="Annotate heatmap values")
    parser.add_argument("--cluster_method", default="average", help="Clustering method")
    parser.add_argument("--colormap", default="viridis", help="Color palette")
    parser.add_argument("--figsize", nargs=2, type=float, default=[12, 10], help="Figure size")
    parser.add_argument("--dpi", type=int, default=300, help="Output DPI")
    parser.add_argument("--vmin", type=float, help="Minimum value for color scale")
    parser.add_argument("--vmax", type=float, help="Maximum value for color scale")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Read input file with multiple encoding attempts
        try:
            df = pd.read_csv(args.input_file, sep=None, engine='python')
        except UnicodeDecodeError:
            df = pd.read_csv(args.input_file, sep=None, engine='python', encoding='latin1')
            
        df.columns = df.columns.str.strip()
        
        ka_ks_col = validate_dataframe(df)
        if not ka_ks_col:
            sys.exit(1)
            
        df = preprocess_data(df, ka_ks_col, args)
        
        # Pivot for heatmaps
        pivot_col = 'Ka_Ks_processed' if 'Ka_Ks_processed' in df else ka_ks_col
        pivoted_df = df.pivot_table(
            index='Gene',
            columns='species_pair',
            values=pivot_col,
            aggfunc='max'
        )
        
        # Create visualizations
        create_heatmap(pivoted_df, args, args.output_dir)
        
        if args.cluster:
            create_clustered_heatmap(pivoted_df, args, args.output_dir)
            
        create_dot_plot(df, ka_ks_col, args, args.output_dir)
        
        # Save processed data
        df.to_csv(f"{args.output_dir}/processed_data.csv", index=False)
        
        print(f"Visualizations saved to {args.output_dir}")
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()