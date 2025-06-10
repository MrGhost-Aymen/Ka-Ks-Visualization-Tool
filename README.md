# Ka/Ks Visualization Tool

A flexible and robust Python script for visualizing Ka/Ks (dN/dS) values across gene pairs from comparative sequence analysis.
its can be used easly with the output of this tool https://github.com/MrGhost-Aymen/Ka-Ks.git
## 📌 Overview

This tool reads a CSV file containing Ka/Ks data and generates publication-quality visualizations including:

- Standard heatmap
- Clustered heatmap (with hierarchical clustering)
- Dot plot showing variation per gene/species pair

The script is designed to handle various column naming conventions, supports log2 transformation, and includes preprocessing of infinite or missing values.

---

## 📦 Requirements

Ensure the following Python packages are installed:

```bash
pip install pandas seaborn matplotlib scipy numpy
```

Optional for enhanced plotting:
- `plotly` (for interactive plots)
- `biopython` (if integrating with sequence analysis pipelines)

---

## 🚀 Usage

```bash
python3 ka_ks_visualizer.py <input_file.csv> [OPTIONS]
```

### 🔧 Example Command

```bash
python3 ka_ks_visualizer.py data.csv \
    --output_dir results \
    --cluster \
    --log_transform \
    --annotate \
    --colormap coolwarm \
    --figsize 16 12 \
    --dpi 600
```

---

## ⚙️ Command-line Options

| Option              | Description |
|---------------------|-------------|
| `--output_dir DIR`  | Output directory (default: `results`) |
| `--cluster`         | Generate clustered heatmap |
| `--log_transform`   | Apply log2 transformation to Ka/Ks values |
| `--annotate`        | Annotate values in heatmaps |
| `--cluster_method METHOD` | Clustering method (e.g., `average`, `ward`) |
| `--colormap NAME`   | Color palette (e.g., `viridis`, `coolwarm`) |
| `--figsize W H`     | Figure size in inches |
| `--dpi N`           | Image resolution in dots-per-inch |
| `--vmin VAL`        | Minimum color scale value |
| `--vmax VAL`        | Maximum color scale value |

---

## 📁 Output Files

All output files are saved in the specified `--output_dir`:

- `heatmap.png`: Main heatmap visualization
- `clustered_heatmap.png`: Hierarchical clustered heatmap *(if enabled)*
- `dot_plot.png`: Stripplot showing gene-level variation
- `processed_data.csv`: Cleaned and processed version of input data

---

## 📄 Input Format

Your input CSV should contain at least these columns:

| Column       | Description |
|--------------|-------------|
| `Gene`       | Gene name or ID |
| `Sequence1`  | First species/strain name |
| `Sequence2`  | Second species/strain name |
| `Ka/Ks`      | Evolutionary divergence ratio (can be named differently; see below) |

### 💡 Supported Ka/Ks Column Names

The script automatically detects one of the following:
- `Ka/Ks`
- `Ka_Ks`
- `KaKs`
- `dN/dS`
- `dn_ds`

---

## 🧪 Example Input (`data.csv`)

```csv
Gene,Sequence1,Sequence2,Ka/Ks,p_value
AT1G01010,Arabidopsis_thaliana,Brassica_rapa,0.45,0.003
AT1G01020,Arabidopsis_thaliana,Oryza_sativa,0.87,0.001
AT1G01010,Brassica_rapa,Oryza_sativa,1.12,0.0001
```

---

## 🧰 Development Tips

- **Add new plots**: Consider adding boxplots, violin plots, or interactive visualizations.
- **Filtering options**: Add CLI flags like `--min_kaks`, `--top_genes`, etc.
- **Logging**: Replace print statements with `logging` module for better control.
- **Testing**: Use `pytest` to test parsing, validation, and plotting logic.

---

## 📬 Feedback & Contributions

Contributions welcome! Please open an issue or submit a PR if you'd like to add features like:

- Interactive plots
- SVG/PDF export
- Multi-sheet Excel output
- Integration with Ka/Ks calculation tools (e.g., PAML, codeml)

---

🌟 **Star this repository** if you find it useful!
