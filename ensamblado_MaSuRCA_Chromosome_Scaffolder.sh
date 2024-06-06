#!/usr/bin/env bash
#SBATCH --mem=200GB
#SBATCH --job-name MaSuRCA_Chromosaome_Scaffolder
#SBATCH --partition=long
#SBATCH --time=1-59:00:20          # total run time limit in HH:MM:SS

#Cargamos los módulos 
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
#Activamos nuestro usuario de mamba para poder acceder a MaSuRCA
mamba init
mamba activate irenegmasurca

cd /home/ireneg/DATOS_DNA-SEQ/procesamiento/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d C*/); do #Recorremos las carpetas de las muestras, en las que se encuentran las lecturas de la secuenciación
        cd $folder #Entramos en la carpeta
        file="${folder:0:18}" #Obtenemos el nombre de la muestra 
        chromosome_scaffolder.sh -r "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome_Assembly_(unmasked)//MagorGY11_1_AssemblyScaffolds.fasta" -q CA/primary.genome.scf.fasta -t 30 -nb #Ejecutamos chromosome scaffolder
        mv MagorGY11_1_AssemblyScaffolds.fasta.primary.genome.scf.fasta.split.reconciled.fa Chromosome_scaffolder #Copiamos el archivo con los scaffodls a una carpeta llamada Chromosome_scaffolder
        cd .. #Salimos de la carpeta de la muestra
done

#-r permite indicar a continuación el genoma de referencia al que se van a referenciar los scaffolds
#-q permite indicar el ensamblaje cuyos scaffolds se van a referenciar al genoma de referencia
#-t indica el número de threads a utilizar
