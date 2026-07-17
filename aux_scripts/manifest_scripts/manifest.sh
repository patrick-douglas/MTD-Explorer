#!/bin/bash

offline_files_folder=

# Entrar na pasta onde serão salvos os arquivos baixados
cd $offline_files_folder

# Arquivo com as URLs a serem baixadas
URLS_FILE="/Kraken2DB_micro/library/bacteria/manifest.list.txt"

# Arquivo para registrar downloads falhados
FAILED_DOWNLOADS="failed_downloads.txt"

# Função para verificar se o arquivo já existe
file_exists() {
    local file=$1
    [ -f "$file" ]
}

# Função para baixar os arquivos com aria2c apenas se ainda não foram baixados
download_files() {
    while read -r url; do
        echo "Verificando: $url"
        filename=$(basename "$url")

        if file_exists "$filename"; then
            echo "Arquivo $filename já existe, pulando download."
        else
            echo "Baixando: $url"
            aria2c --auto-file-renaming=false --continue -x 16 -s 16 -o "$filename" "$url"
            exit_code=$?

            if [ $exit_code -ne 0 ]; then
                echo "Download de $url falhou (Código de saída: $exit_code)"
                echo "$url" >> "$FAILED_DOWNLOADS"
            fi
        fi
    done < "$URLS_FILE"
}

# Executa a função de download
download_files

# Verifica se houve downloads falhados
if [ -s "$FAILED_DOWNLOADS" ]; then
    echo "Downloads que falharam:"
    cat "$FAILED_DOWNLOADS"
else
    echo "Todos os downloads foram concluídos com sucesso!"
fi

