import pandas as pd
from sklearn.cross_decomposition import PLSRegression
import matplotlib.pyplot as plt
import argparse

# Configure argparse to capture command-line arguments
parser = argparse.ArgumentParser(description="Perform PLS-DA on microbiomes and host gene datasets.")
parser.add_argument("-x", "--microbiomes", required=True, help="Path to the Microbiomes.txt file")
parser.add_argument("-y", "--host_gene", required=True, help="Path to the Host_gene.txt file")
parser.add_argument("-o", "--output_pdf", required=True, help="Path to the output PDF file")

args = parser.parse_args()

# Load microbiomes and host gene data
microbiomes_df = pd.read_csv(args.microbiomes, sep='\t', index_col=0)
host_gene_df = pd.read_csv(args.host_gene, sep='\t', index_col=0)

# Fill NaN values with zero to avoid issues in PLS-DA
microbiomes_df = microbiomes_df.fillna(0)
host_gene_df = host_gene_df.fillna(0)

# Ensure samples align across both datasets
common_samples = microbiomes_df.columns.intersection(host_gene_df.columns)
microbiomes_df = microbiomes_df[common_samples].T
host_gene_df = host_gene_df[common_samples].T

# Configure the PLS-DA model with 2 components
pls_da = PLSRegression(n_components=2)

# Fit the model and project data onto PLS components
pls_da.fit(microbiomes_df, host_gene_df)
X_scores, Y_scores = pls_da.transform(microbiomes_df, host_gene_df)

# Plot the results with differentiated sample labels
plt.figure(figsize=(10, 8))
plt.scatter(X_scores[:, 0], X_scores[:, 1], label="Microbiomes", alpha=0.7, c='blue')
plt.scatter(Y_scores[:, 0], Y_scores[:, 1], label="Host Genes", alpha=0.7, c='green')

# Add differentiated labels for each point
for i, sample in enumerate(common_samples):
    plt.text(X_scores[i, 0], X_scores[i, 1], f"{sample}_Microbiome", fontsize=8, color='blue', alpha=0.7)
    plt.text(Y_scores[i, 0], Y_scores[i, 1], f"{sample}_Gene", fontsize=8, color='green', alpha=0.7)

# Labels and title
plt.xlabel("PLS Component 1")
plt.ylabel("PLS Component 2")
plt.title("PLS-DA of Microbiomes and Host Genes with Differentiated Labels")
plt.legend()
plt.tight_layout()

# Save the plot as PDF
plt.savefig(args.output_pdf, format='pdf')
plt.close()
