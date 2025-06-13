#!/bin/bash

# Verifica si se proporcionó un archivo
if [ $# -ne 1 ]; then
    echo "Uso: $0 archivo.fasta"
    exit 1
fi

fasta=$1

# Variables de control
contig=""
length=0

while read -r line; do
    if [[ $line == ">"* ]]; then
        # Si ya hay un contig anterior, imprimimos su longitud
        if [[ -n $contig ]]; then
            echo -e "$contig\t$length"
        fi
        # Guardamos el nuevo nombre del contig y reiniciamos longitud
        contig=${line#>}
        length=0
    else
        # Sumamos la longitud de la línea (nucleótidos)
        line=${line//[[:space:]]/}  # Eliminamos espacios y saltos si los hay
        length=$((length + ${#line}))
    fi
done < "$fasta"

# Imprime el último contig
if [[ -n $contig ]]; then
    echo -e "$contig\t$length"
fi
