#!/usr/bin/env bash
#SBATCH --mem=100GB
#SBATCH --job-name SOAPdenovo2
#SBATCH --partition=long
#SBATCH --time=2-00:00:00          # total run time limit in HH:MM:SS

#Cargamos los módulos y el programa:
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load SOAPdenovo2/r242


cd /home/ireneg/DATOS_DNA-SEQ/originales_descompr/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d */); do  #Recorremos las carpetas de las muestras
        cd $folder #Entramos en la carpeta
        file="${folder:0:18}" #Obtenemos el nombre de la muestra

        #Creamos el archivo de configuración en base a https://github.com/aquaskyline/SOAPdenovo2, adaptándolo a nuestros datos
        echo "#maximal read length" > "SOAPdenovo2/config_file"
        echo "max_rd_len=101" >> "SOAPdenovo2/config_file"
        echo "[LIB]" >> "SOAPdenovo2/config_file"
        echo "#average insert size of the library" >> "SOAPdenovo2/config_file"
        echo "avg_ins=275" >> "SOAPdenovo2/config_file" #Este dato ha sido obtenido tras un análisis de Picard
        echo "#if sequences are forward-reverse of reverse-forward" >> "SOAPdenovo2/config_file"
        echo "reverse_seq=0" >> "SOAPdenovo2/config_file"
        echo "#in which part(s) the reads are used (only contigs, only scaffolds, both contigs and scaffolds, only gap closure)" >> "SOAPdenovo2/config_file"
        echo "asm_flags=3" >> "SOAPdenovo2/config_file"
        echo "#in which order the reads are used while scaffolding" >> "SOAPdenovo2/config_file"
        echo "rank=1" >> "SOAPdenovo2/config_file"
        echo "# cutoff of pair number for a reliable connection (at least 3 for short insert size)" >> "SOAPdenovo2/config_file"
        echo "pair_num_cutoff=3" >> "SOAPdenovo2/config_file"
        echo "#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)" >> "SOAPdenovo2/config_file"
        echo "map_len=32" >> "SOAPdenovo2/config_file"
        echo "#paired-end fastq files, read 1 file should always be followed by read 2 file" >> "SOAPdenovo2/config_file"
        
        echo q1=out_fastp_${file}_1_paired.fq >> "SOAPdenovo2/config_file" 
        echo q2=out_fastp_${file}_2_paired.fq >> "SOAPdenovo2/config_file"
        
        echo q=out_fastp_${file}_1_unpaired.fq >> "SOAPdenovo2/config_file"
        echo q=out_fastp_${file}_2_unpaired.fq >> "SOAPdenovo2/config_file"
                
        SOAPdenovo-63mer all -s SOAPdenovo2/config_file -o SOAPdenovo2/K21/output -K 21 #Ejecutamos SOAPdenovo con un tamaño de kmer de 21
        SOAPdenovo-63mer all -s SOAPdenovo2/config_file -o SOAPdenovo2/K33/output -K 33 #Ejecutamos SOAPdenovo con un tamaño de kmer de 33
        SOAPdenovo-63mer all -s SOAPdenovo2/config_file -o SOAPdenovo2/K55/output -K 55 #Ejecutamos SOAPdenovo con un tamaño de kmer de 55
        
        cd ..
done

#En el comando de ejecución de SOAPdenovo2:
#all indica que se deben ejecutar todas las etapas del ensamblaje
#-s indica la localización del archivo de configuración, en el que están especificados los parámetros a utilizar
#-o indica dónde se debe guardadr el output
#-k indica el tamaño de kmer utilizado para construir el grafo de Bruijn. En nuestor caso, se ha ejecutado el programa con 3 tamaños para optimizar el proceso. 
