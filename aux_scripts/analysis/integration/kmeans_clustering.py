import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import argparse
from sklearn.preprocessing import StandardScaler

# Configuração do argparse para capturar argumentos da linha de comando
parser = argparse.ArgumentParser(description="Perform K-means clustering on microbiomes and host gene datasets.")
parser.add_argument("-x", "--microbiomes", required=True, help="Path to the Microbiomes.txt file")
parser.add_argument("-y", "--host_gene", required=True, help="Path to the Host_gene.txt file")
parser.add_argument("-o", "--output_pdf", required=True, help="Path to the output PDF file")
parser.add_argument("-k", "--clusters", type=int, default=3, help="Number of clusters for K-means")

args = parser.parse_args()

# Carregar os dados de microbiomas e genes hospedeiros
microbiomes_df = pd.read_csv(args.microbiomes, sep='\t', index_col=0)
host_gene_df = pd.read_csv(args.host_gene, sep='\t', index_col=0)

# Preencher valores NaN com zero para evitar problemas
microbiomes_df = microbiomes_df.fillna(0)
host_gene_df = host_gene_df.fillna(0)

# Garantir que as amostras estejam alinhadas
common_samples = microbiomes_df.columns.intersection(host_gene_df.columns)
data_combined = pd.concat([microbiomes_df[common_samples].T, host_gene_df[common_samples].T], axis=1)

# Normalizar os dados
scaler = StandardScaler()
data_normalized = scaler.fit_transform(data_combined)

# Aplicar K-means clustering
kmeans = KMeans(n_clusters=args.clusters, random_state=42)
kmeans.fit(data_normalized)
labels = kmeans.labels_

# Plotar os clusters
plt.figure(figsize=(10, 8))
plt.scatter(data_normalized[:, 0], data_normalized[:, 1], c=labels, cmap='viridis', alpha=0.7)
plt.xlabel("Feature 1 (Normalized)")
plt.ylabel("Feature 2 (Normalized)")
plt.title(f"K-means Clustering of Microbiomes and Host Genes (k={args.clusters})")
plt.colorbar(label='Cluster Label')
plt.tight_layout()

# Salvar o gráfico como PDF
plt.savefig(args.output_pdf, format='pdf')
plt.close()
