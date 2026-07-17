import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Configure argparse to capture command-line arguments
parser = argparse.ArgumentParser(description="Generate a heatmap of mean abundances/expressions for microbiomes and host genes.")
parser.add_argument("-x", "--microbiomes", required=True, help="Path to the Microbiomes.txt file")
parser.add_argument("-y", "--host_gene", required=True, help="Path to the Host_gene.txt file")
parser.add_argument("-o", "--output_pdf", required=True, help="Path to the output PDF file")
parser.add_argument("--x_dataset_label", default="Microbiomes", help="Label for the x dataset (Microbiomes)")
parser.add_argument("--y_dataset_label", default="Host Gene", help="Label for the y dataset (Host Gene)")

args = parser.parse_args()

# Load microbiomes and host gene data
microbiomes_df = pd.read_csv(args.microbiomes, sep='\t', index_col=0).fillna(0)
host_gene_df = pd.read_csv(args.host_gene, sep='\t', index_col=0).fillna(0)

# Calculate the mean of each sample for microbiomes and host genes
microbiomes_mean = microbiomes_df.mean(axis=0)
host_gene_mean = host_gene_df.mean(axis=0)

# Create a DataFrame with the means for visualization
mean_df = pd.DataFrame({f'{args.x_dataset_label} Mean': microbiomes_mean, f'{args.y_dataset_label} Mean': host_gene_mean})

# Plot the heatmap of means
plt.figure(figsize=(10, 8))
sns.heatmap(mean_df.T, cmap="viridis", annot=True, fmt=".2f", cbar_kws={'label': 'Mean Expression/Abundance'})

# Labels
plt.xlabel("Samples")
plt.ylabel("Mean per Group")
plt.title("Heatmap of Means for Microbiomes and Host Genes")
plt.tight_layout()

# Save the plot as a PDF
plt.savefig(args.output_pdf, format='pdf')
