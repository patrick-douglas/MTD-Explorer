#!/bin/bash
# Check R package versions in a nice table

# Cores
w=$(tput sgr0)
g=$(tput setaf 2)
r=$(tput setaf 1)

# Versão do R
R_ver=$(R --version | grep version | grep R | awk '{print $3}')

# Pacotes a checar
pkgs=( lattice MASS mnormt nlme GPArotation psych foreign R.methodsS3 R.oo rtf psychTools XICOR mclust BiocManager preprocessCore remotes EnvStats Hmisc eva Matrix )

# Cabeçalho da tabela
echo "${g}╔═════════════════╦═══════════════╗"
echo "║R                ║ ${w}$R_ver${g}         ║"
echo "║Conda Environment║ ${w}halla0820${g}     ║"
echo "╠═════════════════╬═══════════════╣"

# Loop pelos pacotes
for pkg in "${pkgs[@]}"; do
    ver=$(R --no-restore -e "packageVersion(\"${pkg}\")" 2>/dev/null | \
           grep '\[1\]' | awk '{print $2}' | sed -r 's/^.{1}//; s/.$//')
    
    # Se não estiver instalado
    if [ -z "$ver" ]; then
        ver="${r}not installed${g}"
    fi

    printf "║ %-15s ║ ${w}%-13s${g} ║\n" "$pkg" "$ver"
done

# Rodapé da tabela
echo "╚═════════════════╩═══════════════╝"
echo "${w}"

