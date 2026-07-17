# Definir um espelho CRAN padrão (pode ajustar para um mirror mais próximo, se necessário)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Função para verificar e instalar pacotes se necessário
install_if_missing <- function(package) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package, dependencies = TRUE)
        library(package, character.only = TRUE)
    }
}

# Verificar e instalar pacotes necessários
install_if_missing("VennDiagram")
install_if_missing("argparser")

# Função para extrair a lista de espécies dos arquivos .krona
extract_species_list <- function(krona_files) {
    combined_species <- c()
    for (krona_file in krona_files) {
        temp_file <- paste0(krona_file, "_species_list.txt")
        grep_command <- sprintf("grep 's__' %s | awk -F 's__' '{print $2}' | sort | uniq > %s", krona_file, temp_file)
        system(grep_command)
        species <- scan(temp_file, what = "", sep = "\n", quiet = TRUE)
        combined_species <- unique(c(combined_species, species))
        file.remove(temp_file)
    }
    return(combined_species)
}

# Criar o parser de argumentos
p <- arg_parser("Generate a Venn diagram based on multiple input krona files.")
p <- add_argument(p, "--krona_files1", help="Paths to first group krona file(s)", nargs = Inf)
p <- add_argument(p, "--krona_files2", help="Paths to second group krona file(s)", nargs = Inf)
p <- add_argument(p, "--group_label1", help="Label for first group")
p <- add_argument(p, "--group_label2", help="Label for second group")

# Parse argumentos
args <- parse_args(p)

# Extrair espécies dos dois grupos de arquivos krona
species_group1 <- extract_species_list(args$krona_files1)
species_group2 <- extract_species_list(args$krona_files2)

# Calcular espécies únicas e compartilhadas
unique_species_group1 <- setdiff(species_group1, species_group2)       # Espécies únicas para o grupo 1
unique_species_group2 <- setdiff(species_group2, species_group1)       # Espécies únicas para o grupo 2
shared_species <- intersect(species_group1, species_group2)            # Espécies compartilhadas

# Calcular os totais corrigidos para cada grupo
group1_unique <- length(unique_species_group1)          # Espécies exclusivas para o grupo 1
group2_unique <- length(unique_species_group2)          # Espécies exclusivas para o grupo 2
intersection <- length(shared_species)                  # Espécies compartilhadas

# Definir os raios dos círculos com base em fatores ajustados
radius1 <- 0.15 + (group1_unique + intersection) / 2000
radius2 <- 0.15 + (group2_unique + intersection) / 2000

# Calcular as coordenadas automáticas para posicionamento
# Círculo 1 (esquerda)
center_x1 <- 0.4
center_y1 <- 0.55

# Círculo 2 (direita)
center_x2 <- 0.6
center_y2 <- 0.55

# Interseção (centro)
center_intersection_x <- (center_x1 + center_x2) / 2
center_intersection_y <- center_y1  # Eles compartilham o mesmo eixo Y

# Definir deslocamento adicional horizontal
offset <- 0.05  # Pequeno deslocamento para afastar um pouco mais para as extremidades

# Desenhar os círculos manualmente com tamanhos ajustados
draw_venn_diagram <- function() {
    # Desenhar o círculo do grupo 1 (esquerda)
    grid.circle(x = center_x1, y = center_y1, r = radius1, gp = gpar(fill = "#fc9272", alpha = 0.5, col = "black"))
    
    # Desenhar o círculo do grupo 2 (direita)
    grid.circle(x = center_x2, y = center_y2, r = radius2, gp = gpar(fill = "#a1d99b", alpha = 0.5, col = "black"))
    
    # Centralizar os números nas áreas correspondentes com deslocamento horizontal
    grid.text(label = as.character(group1_unique), x = center_x1 - radius1 / 2 - offset, y = center_y1, gp = gpar(fontsize = 18, fontface = "plain"))  # Posição no vermelho
    grid.text(label = as.character(group2_unique), x = center_x2 + radius2 / 2 + offset, y = center_y2, gp = gpar(fontsize = 18, fontface = "plain"))  # Posição no verde
    grid.text(label = as.character(intersection), x = center_intersection_x, y = center_intersection_y, gp = gpar(fontsize = 18, fontface = "plain"))   # Posição no centro (amarelo)
    
    # Título no topo (maior que os números e rótulos)
    grid.text("Total number of species detected in the microbiome", y = 0.98, gp = gpar(fontsize = 22, fontface = "bold"))
    
    # Rótulos dos grupos abaixo dos círculos (maiores que os números, mas menores que o título)
    grid.text(args$group_label1, x = 0.32, y = 0.15, gp = gpar(fontsize = 20, fontface = "bold"))
    grid.text(args$group_label2, x = 0.68, y = 0.15, gp = gpar(fontsize = 20, fontface = "bold"))
}

# Salvar o diagrama em vários formatos, com o canvas aumentado
output_base <- sprintf("venn_diagram_%s_vs_%s", args$group_label1, args$group_label2)

# Aqui está o ajuste para o PNG, com um canvas maior para evitar corte do título
png(paste0(output_base, ".png"), width = 1800, height = 2000, res = 300)  # Aumentar a altura do canvas
draw_venn_diagram()
dev.off()

pdf(paste0(output_base, ".pdf"), width = 10, height = 10)  # Canvas grande para PDF
draw_venn_diagram()
dev.off()

svg(paste0(output_base, ".svg"), width = 10, height = 10)  # Canvas grande para SVG
draw_venn_diagram()
dev.off()

# Remover o arquivo Rplots.pdf se ele existir
if (file.exists("Rplots.pdf")) {
    file.remove("Rplots.pdf")
}