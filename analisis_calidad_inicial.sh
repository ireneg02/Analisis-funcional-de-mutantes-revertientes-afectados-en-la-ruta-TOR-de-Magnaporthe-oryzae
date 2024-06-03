#!/usr/bin/env bash
#SBATCH --mem=32GB
#SBATCH --job-name analisis_calidad_inicial
#SBATCH --partition=fast
#SBATCH --time=02:00:00

#Importamos los módulos del clúster y el programa FastQC
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load FastQC/0.11.8-Java-1.8

#
fastqc -o ~/DATOS_DNA-SEQ/dnaseq/results/resultados_originales -t 6 *.fastq.gz

#Activamos nuestro usuario de mamba para poder acceder a MultiQC
mamba init
mamba activate irenegmasurca

#Ejecutamos MultiQC
multiqc /home/ireneg/DATOS_DNA-SEQ/FASTQC/results/resultados_originales/ -o multiqc_fastqc_inicial -p
#-o indica el nombre del directorio en el que se va a guardar el output
#-p se usa para exportar los gráficos 