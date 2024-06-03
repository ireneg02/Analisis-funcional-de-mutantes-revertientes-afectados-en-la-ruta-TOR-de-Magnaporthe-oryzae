#!/usr/bin/env bash
#SBATCH --mem=100GB
#SBATCH --job-name filtrado_trimmomatic
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00

#Cargamos los modulos y el programa:
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load Java/1.8.0_212
module load Trimmomatic/0.38-Java-1.8

cd /home/ireneg/DATOS_DNA-SEQ/originales #Nos movemos a la carpeta con los archivos de la secuenciaci�n

#Los archivos sobre los que vamos a trabajar est�n en formato fasta.gz. Tenemos un archivo con las lecturas forward que termina en 1 y otro con las lecturas reverse que termina en 2. 
#As�, todos los archivos siguen la estructura "numero_muestra_1.fasta.gz" para el archivo con lecturas forward y "numero_muestra_2.fasta.gz" para el archivo con las lecturas reverse.

for file_1 in $(ls *1.fastq.gz); do #Dentro de la carpeta en la que est�n todos los archivos de la secuenciaci�n, recorremos todos los archivos forward (terminan en 1.fastq.gz)
        file_2="${file_1:0:19}2.fastq.gz" #Obtenemos el nombre del archivo reverse que va en pareja con el que estamos procesando en el bucle
        file="${file_1:0:19}" #Obtenemos el nombre del archivo de la muestra (eliminando el 1 o 2 del final y la extensi�n)
        java -jar /home/sbodi/miniconda3/envs/SqueezeMeta/SqueezeMeta/bin/trimmomatic-0.38.jar PE -phred33 -trimlog trimLog_${file} -summary statsSummary_${file} $file_1 $file_2 output_${file_1}_paired.fq.gz output_${file_1}_unpaired.fq.gz output_${file_2}_paired.fq.gz output_${file_2}_unpaired.fq.gz ILLUMINACLIP:mis_adaptadores.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 #Ejecutamos Trimmomatic (especificaci�n de los comandos a continuaci�n)
done


#PE significa que vamos a realizar un an�lisis que es un Paired Ends

#-phred33 indica la escala de calidad de las secuencia por la que vamos a filtrar

#-trimLog permite crear un archivo log con informaci�n sobre las lecturas procesadas que contiene: 
#el nombre de la lectura
#la longitud de la secuencia superviviente
#la posici�n de la primera base superviviente, es decir, la cantidad recortada desde el principio
#la posici�n de la �ltima base superviviente en la lectura original
#la cantidad recortada del final

#-summary especifica que se debe generar un archivo con informaci�n sobre el proceso de filtrado. Este archivo incluye el n�mero de pares de lectura de entrada, de lecturas supervivientes (y el porcentaje), de lecturas eliminadas (y el porcentaje).

#A continuacion tenemos el archivo forward a procesar y luego el reverse. Despu�s se especifican los archivos en los que se van a almacenar las lecturas forward que est�n apareadas, las forward desapareadas, las reverse apareadas y las reverse desapareadas en ese orden.

#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<keepBothReads>
#fastaWithAdaptersEtc: especifica la ruta a un archivo fasta que contiene todos los adaptadores.
#seedMismatches: especifica el n�mero m�ximo de desajustes que permitir�n realizar una comparaci�n completa.
#palindromeClipThreshold: especifica la precisi�n que debe tener la coincidencia entre las dos lecturas �ligadas al adaptador� para la alineaci�n de lecturas palindr�micas PE.
#simpleClipThreshold: especifica la precisi�n de la coincidencia entre cualquier secuencia de adaptador con una lectura. 
#keepBothReads: Cuando toma el valor "True", indica que un par de lecturas emparejadas deben mantenerse incluso si solo una de ellas tiene un adaptador recortado.
#ILLUMINACLIP:mis_adaptadores.fa:2:30:10:2:True indica que se tomen los adaptadores del archivo "mis_adaptadores", se permiten hasta 2 desajustes en la comparaci�n de las lecturas y los adaptadores, requerimos una coincidencia de al menos 30 bases para considerar y recortar las lecturas palindr�micas, necesitamos una coincidencia de al menos 10 bases para considerar y recortar el adaptador de la lectura, cualquier secuencia de adaptador que coincida en al menos 2 bases ser� recortada y que un par de lecturas emparejadas deben mantenerse incluso si solo una de ellas tiene un adaptador recortado

#LEADING corta las bases al inicio de una lectura cuadno estas est�n por debajo de una calidad umbral.
#TRAILING hace lo mismo que LEADING pero aplic�ndose al final de las lecturas.
#Por lo tanto, LEADING:3 y TRAILING:3 hacen que se eliminen las bases al principio y al final de una lectura que tienen una calidad menor a 3

#SLIDINGWINDOW realiza un recorte de ventana deslizante, cortando una vez que la calidad media dentro de la ventana cae por debajo de un umbral. En particular SLIDINGWINDOW:4:20 comprueba que la calidad promedio en una ventana de 4 bases no sea menor de 20. El est�ndar suele ser usar una calidad que no caiga de 15, pero como nuestros datos tienen una alta calidad, hemos podido subir el umbral para aplicar un filtrado m�s exhaustivo. 

#MINLEN se usa para especificar la longitud m�nima que deben tener las lecturas para mantenerse. Con MINLEN:36 se mantienen solo las secuencias con una longitud mayor a 36 pares de bases. Este valor es el usado como est�ndar. 

#Documentaci�n en https://github.com/usadellab/Trimmomatic