#!/usr/bin/env bash
#SBATCH --mem=100GB
#SBATCH --job-name coverage_SAMBEDtools
#SBATCH --partition=medium
#SBATCH --time=24:00:00          # total run time limit in HH:MM:SS

#Cargamos los módulos y programas:
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load BEDTools/2.28.0-GCC-8.2.0-2.31.1
module load SAMtools/1.9-GCC-8.2.0-2.31.1

cd /home/ireneg/DATOS_DNA-SEQ/procesamiento/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d C*/); do #Recorremos las carpetas de las muestras
        cd $folder #Entramos en la carpeta de la muestra
        file="${folder:0:18}" #Obtenemos el nombre de la muestra
        samtools depth -a BWA/sorted_bwa_map_${file}_merge_PE_SE_group.bam > BWA/coverage_${file}.txt #Calculamos la profundidad de cobertura del archivo BAM y guardamos la salida en un archivo de texto
        awk 'BEGIN {OFS="\t"} $3 == 0 {print $1, $2, $2 + 1}' BWA/coverage_${file}.txt > BWA/zero_coverage_${file}.bed #Filtramos las posiciones con cobertura cero y las guardarmos en un archivo BED
        bedtools merge -i BWA/zero_coverage_${file}.bed > BWA/merged_zero_coverage_${file}.bed #Fusionamos las regiones contiguas con cobertura cero usando bedtools
        awk 'BEGIN {OFS="\t"} {print $1, "bedtools", "zero_coverage", $2, $3, ".", ".", ".", "coverage=0"}' BWA/merged_zero_coverage_${file}.bed > BWA/zero_coverage_${file}.gff #Convertimos las regiones fusionadas a formato GFF
        awk '{if ($5 - $4 >= 100) print $0}' BWA/zero_coverage_${file}.gff > BWA/zero_coverage_l100_${file}.gff #Filtramos las regiones GFF para que solo queden las de al menos 100 bases de longitud.
        cd .. #Salimos de la carpeta de la muestra
done

#En el comando de samtools depth:
#-a especifica que se incluyan todas las posiciones en el archivo de salida, incluidas las de profundidad de coverage 0. 
#A continuación, se indica el nombre del archivo BAM (derivado del mapeado) a procesar
#Finalmente, tras >, se indica el archivo en el que se va a almacenar la salida

#En uso de awk para filtrar las posiciones de profundidad cero, se seleccionan las filas cuya tercera columna (en la que se almacena la profundidad de coverage) es 0. El resultado se almacena en un archivo bed

#BEDtools merge se usa para fusionar intervalos en el archivo BED generado anteriormente. 
#-i indica el archivo de entrada
#Después de ">", se indica el archivo generado con las regiones de profundidad cero fusionadas

#La segunda vez que empleamos awk es para reorganizar las columnas para adaptarlas al formato gff, generando como output un archivo gff 

#Finalmente, el último comando de awk se emplea para filtrar las filas donde la diferencia entre la quinta (que almacena la posición final de la región con profundidad cero) y la cuarta columna (que almacena la posición inicial) es mayor o igual a 100. Así, estamos eliminando todas las regiones que tengan un tamaño menor a 100pb. 
