#!/usr/bin/env bash
#SBATCH --mem=32GB
#SBATCH --job-name analisis_calidad_inicial
#SBATCH --partition=fast
#SBATCH --time=02:00:00

#Importamos los módulos del clúster y el programa FastQC
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load FastQC/0.11.8-Java-1.8

#Importamos los módulos del clúster y el programa FastQC
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load FastQC/0.11.8-Java-1.8

#Ejecutamos el análisis de FastQC
fastqc -o ~/DATOS_DNA-SEQ/originales/FASTQC -t 6 *.fastq.gz
#-t especifica el número de threads que se deben usar
#A continuación, se especifica el/los archivos sobre los que se debe hacer el análisis. En este caso, usando wildcards podemos agilizar que se realice en todos los archivos de la secuenciación, sin tener que especificar todos los nombres. 

#Activamos nuestro usuario de mamba para poder acceder a MultiQC
mamba init
mamba activate irenegmasurca

#Ejecutamos MultiQC
multiqc /home/ireneg/DATOS_DNA-SEQ/originales/FASTQC -o multiqc_fastqc_inicial -p
#-o indica el nombre del directorio en el que se va a guardar el output
#-p se usa para exportar los gráficos 
