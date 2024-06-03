#!/usr/bin/env bash
#SBATCH --mem=100GB
#SBATCH --job-name filtrado_fastp
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00          # total run time limit in HH:MM:SS

#Cargamos los módulos y el programa:
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
PATH="$PATH:/home/ireneg/"

cd /home/ireneg/DATOS_DNA-SEQ/procesamiento/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d */); do #Recorremos las carpetas de las muestras, en las que se encuentran las lecturas de la secuenciación
        cd $folder #Entramos en la carpeta
        file="${folder:0:18}" #Obtenemos el nombre de la muestra
        file_1="${file}_1.fastq.gz" #Obtenemos el nombre del archivo con las lecturas forward
        file_2="${file}_2.fastq.gz" #Obtenemos el nombre del archivo con las lecturas reverse  
        fastp -i $file_1 -I $file_2 --out1 out_fastp_${file}_1_paired.fq.gz --out2 out_fastp_${file}_2_paired.fq.gz --unpaired1 out_fastp_${file}_1_unpaired.fq.gz --unpaired2 out_fastp_${file}_2_unpaired.fq.gz --failed_out out_fastp_${file}_failed -q 15 -l 36 --detect_adapter_for_pe --cut_right -M 20 -u 40 -j ${file}.fastp.json -h ${file}.fastp.html #Ejecutamos el filtrado con Fastp. Comandos explicados a continuación.
        
        #Descomprimimos los archivos manteniendo los originales
        gunzip -k out_fastp_${file}_1_paired.fq.gz
        gunzip -k out_fastp_${file}_2_paired.fq.gz
        gunzip -k out_fastp_${file}_1_unpaired.fq.gz
        gunzip -k out_fastp_${file}_2_unpaired.fq.gz
        cd ..
done


# -i permite introducir el file con datos forward y -I los del reverse
#--out1 indica el nombre del archivo en el que se almacenan las secuencias filtradas forward emparejadas, --out2 de las filtradas reverse emparejadas, --unpaired1 de las secuencias filtradas forward desemparejadas, --unpaired2 filtradas reverse desemparejadas y --failed_out las secuencias que no han superado el filtrado
#-q indica el nivel mínimo de calidad de las bases individuales, por lo que las bases con menor calidad se eliminan
#-l establece la longitud mínima de las secuencias
#--detect_adapter_for_pe sirve para que identifique los adaptadores automáticamente y los filtre
#--cut_right actúa como una ventana deslizante que mide la calidad media de un número de residuos especificados (el default son 4 residuos) y, si es menor que un umbral especificado, la elimina. El umbral de calidad por defecto es 20, pero se puede modificar con el comando -M. Nosotros lo mantenemos, pues es el mismo criterio que habíamos utilizado con Trimmomatic. 
#-u indica el porcentaje de N (bases indefinidas) que se permiten. En este caso solo se permite un 40% de las N
#-j  permite definir el nombre del archivo json que se va a generar
#-h  permite definir el nombre del html que se va a generar
