#!/usr/bin/env bash
#SBATCH --mem=200GB
#SBATCH --job-name BUSCO
#SBATCH --partition=long
#SBATCH --time=3-0:00:00          # total run time limit in HH:MM:SS

#Cargamos los módulos y el programa:
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load BUSCO/5.1.2-foss-2019a

cd /home/ireneg/DATOS_DNA-SEQ/procesamiento/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d C*/); do #Recorremos las carpetas de las muestras, en las que se encuentran las lecturas de la secuenciación
        cd $folder #Entramos en la carpeta
        busco -i spades_output/scaffolds.fasta --auto-lineage -o BUSCO_SPAdes -m genome #ejecutamos el análisis BUSCO en los scaffolds generados por SPAdes
        busco -i SOAPdenovo2/K33/output.scafSeq --auto-lineage -o BUSCO_SOAPdenovo2 -m genome #ejecutamos el análisis BUSCO en los scaffolds generados por SOAPdenovo2
        busco -i CA/primary.genome.scf.fasta --auto-lineage -o BUSCO_MaSuRCA -m genome #ejecutamos el análisis BUSCO en los scaffolds generados por MASuRCA
        busco -i Chromosome_scaffolder/MagorGY11_1_AssemblyScaffolds.fasta.primary.genome.scf.fasta.split.reconciled.fa --auto-lineage -o BUSCO_Chrom_scaff -m genome #ejecutamos el análisis BUSCO en los scaffolds generados por Chromosome scaffolder a partir del output de MaSuRCA
        cd ..
done

#En la ejecución de BUSCO:
#-i especifica el archivo de entrada resultante del ensamblado, que contiene los scaffolds
#--auto-lineage indica que se debe seleccionar automáticamente el linaje más adecuado en base al contenido del archivo de entrada
#-o especifica el directorio de salida donde se guardarán los resultados del análisis
#-m indica el modo de análisis. Si a continuación indicamos "genome", significa que se va a evaluar un ensamblaje de genoma completo.
