#!/usr/bin/env bash
#SBATCH --mem=200GB
#SBATCH --job-name bwa
#SBATCH --partition=long
#SBATCH --time=4-00:00:00

#Cargamos los módulos y el programa:
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load BWA/0.7.17-GCC-8.2.0-2.31.1
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load picard/2.26.10-Java-15

#En primer lugar, ejecutamos el programa mapeando las lecturas al genoma de referencia Guy11
cd "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)" #Nos movemos a la carpeta que contiene el genoma de referencia de Guy11
bwa index MagorGY11_1_AssemblyScaffolds.fasta #Creamos el índice del genoma de referencia Guy11

#Se mapean las lecturas originales de la secuenciación
cd /home/ireneg/DATOS_DNA-SEQ/originales_descompr/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d C*/); do #Recorremos las carpetas de las muestras, en las que se encuentran las lecturas originales sin comprimir
        cd $folder #Entramos en la carpeta
        file="${folder:0:18}" #Obtenemos el nombre de la muestra
        file_1="${file}_1.fastq" #Obtenemos el nombre del archivo con las lecturas forward
        file_2="${file}_2.fastq" #Obtenemos el nombre del archivo con las lecturas reverse     
        bwa mem "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" ${file_1} ${file_2}  > BWA/bwa_map_${file}_unfiltered_Guy11.sam #Ejecutamos el mapeado de las lecturas originales. Los resultados se almacenan en una carpeta llamada "BWA"        
        cd ..#Salimos de la carpeta de la muestra
done

#Se repite el proceso para mapear las lecturas resultantes del filtrado con Fastp
cd /home/ireneg/DATOS_DNA-SEQ/procesamiento/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d */); do #Recorremos las carpetas de las muestras, en las que se encuentran las lecturas emparejadas y sin emparejar resultantes del filtrado con Fastp
        cd $folder #Entramos en la carpeta
        file="${folder:0:18}" #Obtenemos el nombre de la muestra
        file_1="${file}_1.fastq.gz" #Obtenemos el nombre del archivo con las lecturas forward
        file_2="${file}_2.fastq.gz" #Obtenemos el nombre del archivo con las lecturas reverse  
        bwa mem "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" out_fastp_${file}_1_paired.fq out_fastp_${file}_2_paired.fq> BWA/bwa_map_${file}_PE_Guy11.sam #Se mapean las lecturas emparejadas
        bwa mem "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" out_fastp_${file}_1_unpaired.fq > BWA/bwa_map_${file}_SE_1_Guy11.sam #Se mapean las lecturas forward sin pareja
        bwa mem "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" out_fastp_${file}_2_unpaired.fq > BWA/bwa_map_${file}_SE_2_Guy11.sam #Se mapean las lecturas reverse sin pareja
        cd ..#Salimos de la carpeta de la muestra
done
#Al filtrar las lecturas con fastp, algunas de las lecturas quedan desemparejadas, por lo que se deben mapear al genoma de manera independiente

              
        


#En segundo lugar, ejecutamos el programa mapeando las lecturas al genoma de referencia 70-15
cd "/home/ireneg/DATOS_DNA-SEQ/originales/ncbi_dataset/data/GCF_000002495.2/" #Nos movemos a la carpeta que contiene el genoma de referencia de 70-15
bwa index GCF_000002495.2_MG8_genomic.fna #Construimos el índice del genoma de referencia 70-15

#Se mapean las lecturas originales de la secuenciación
cd /home/ireneg/DATOS_DNA-SEQ/originales_descompr/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d C*/); do #Recorremos las carpetas de las muestras, en las que se encuentran las lecturas originales sin comprimir
        cd $folder #Entramos en la carpeta
        file="${folder:0:18}" #Obtenemos el nombre de la muestra
        file_1="${file}_1.fastq" #Obtenemos el nombre del archivo con las lecturas forward
        file_2="${file}_2.fastq" #Obtenemos el nombre del archivo con las lecturas reverse     
        bwa mem "/home/ireneg/DATOS_DNA-SEQ/originales/ncbi_dataset/data/GCF_000002495.2/GCF_000002495.2_MG8_genomic.fna" ${file_1} ${file_2}  > BWA/bwa_map_${file}_unfiltered_70-15.sam #Ejecutamos el mapeado de las lecturas originales. Los resultados se almacenan en una carpeta llamada "BWA"        
        cd ..#Salimos de la carpeta de la muestra
done

#Se repite el proceso para mapear las lecturas resultantes del filtrado con Fastp
cd /home/ireneg/DATOS_DNA-SEQ/procesamiento/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d */); do #Recorremos las carpetas de las muestras, en las que se encuentran las lecturas emparejadas y sin emparejar resultantes del filtrado con Fastp
        cd $folder #Entramos en la carpeta
        file="${folder:0:18}" #Obtenemos el nombre de la muestra
        file_1="${file}_1.fastq.gz" #Obtenemos el nombre del archivo con las lecturas forward
        file_2="${file}_2.fastq.gz" #Obtenemos el nombre del archivo con las lecturas reverse  
        bwa mem "/home/ireneg/DATOS_DNA-SEQ/originales/ncbi_dataset/data/GCF_000002495.2/GCF_000002495.2_MG8_genomic.fna" out_fastp_${file}_1_paired.fq out_fastp_${file}_2_paired.fq> BWA/bwa_map_${file}_PE_70-15.sam #Se mapean las lecturas emparejadas
        bwa mem "/home/ireneg/DATOS_DNA-SEQ/originales/ncbi_dataset/data/GCF_000002495.2/GCF_000002495.2_MG8_genomic.fna" out_fastp_${file}_1_unpaired.fq > BWA/bwa_map_${file}_SE_1_70-15.sam #Se mapean las lecturas forward sin pareja
        bwa mem "/home/ireneg/DATOS_DNA-SEQ/originales/ncbi_dataset/data/GCF_000002495.2/GCF_000002495.2_MG8_genomic.fna" out_fastp_${file}_2_unpaired.fq > BWA/bwa_map_${file}_SE_2_70-15.sam #Se mapean las lecturas reverse sin pareja
        cd ..#Salimos de la carpeta de la muestra
done