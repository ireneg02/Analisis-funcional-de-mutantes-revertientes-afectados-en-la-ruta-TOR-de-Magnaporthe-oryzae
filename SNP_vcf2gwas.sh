#!/usr/bin/env bash
#SBATCH --mem=100GB
#SBATCH --job-name vcf2gwas
#SBATCH --partition=medium
#SBATCH --time=5:00:00          # total run time limit in HH:MM:SS

#Cargamos los módulos y el programa:
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load BCFtools/1.9-GCC-8.2.0-2.31.1

cd /home/ireneg/DATOS_DNA-SEQ/procesamiento/gwas/ #Nos desplazamos a la carpeta que contiene los archivos con los SNP tras el filtrado de cada muestra

#Convertimos el formato para poder trabajar con bcftools. Asociamos el número de muestra con su nombre. 
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_12nfDI.vcf> Guy11_0.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_24nfDI.vcf> 2D4_0.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_36nfDI.vcf> ID4.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_48nfDI.vcf> ID39.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_59nfDI.vcf> Guy11_1.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_60nfDI.vcf> ID42.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_71nfDI.vcf> Guy11_2.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_72nfDI.vcf> ID46A.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_83nfDI.vcf> 2D4_1.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_84nfDI.vcf> ID46B.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_95nfDI.vcf> 2D4_2.vcf.gz
bgzip -c SNP_qualfilter_BWA_CBMUUANXX_6_96nfDI.vcf> ID49.vcf.gz

#Indexamos los archivos
bcftools index Guy11_0.vcf.gz
bcftools index 2D4_0.vcf.gz
bcftools index ID4.vcf.gz
bcftools index ID39.vcf.gz
bcftools index Guy11_1.vcf.gz
bcftools index ID42.vcf.gz
bcftools index Guy11_2.vcf.gz
bcftools index ID46A.vcf.gz
bcftools index 2D4_1.vcf.gz
bcftools index ID46B.vcf.gz
bcftools index 2D4_2.vcf.gz
bcftools index ID49.vcf.gz

#Fusionamos los vcf de todas las muestras en un único archivo
bcftools merge Guy11_0.vcf.gz 2D4_0.vcf.gz ID4.vcf.gz ID39.vcf.gz Guy11_1.vcf.gz ID42.vcf.gz Guy11_2.vcf.gz ID46A.vcf.gz 2D4_1.vcf.gz ID46B.vcf.gz 2D4_2.vcf.gz ID49.vcf.gz --force-samples -o merged_SNP.vcf.gz 
#--force-samples fuerza la fusión incluso si hay conflictos en los nombres de las muestras
#-o indica el nombre del archivo output

bcftools sort merged_SNP.vcf.gz -O v -o sorted_merged_SNP.vcf
#-O v: especifica que la salida debe ser en formato VCF sin comprimir
#-o indica el nombre del archivo de salida

awk '
BEGIN { FS=OFS="\t" }
{
    if ($0 !~ /^#/) {
        for (i=10; i<=NF; i++) {
            if ($i == "./.:.:.:.:.") {
                $i = "0/0"; #Se sustituyen los campos "./.:.:.:.:.", característicos de que una variante no esté presente en la muestra, por 0/0, que es el genotipo usado para indicar que no se observan diferencias respecto del genoma de referencia
            } else {
                split($i, geno_info, ":");
                $i = geno_info[1];  # Conservamos solo la información del genotipo (GT)
            }
        }
        $8 = ".";  # Reemplazamos el contenido de la octava columna por "."
        $9 = "GT"; # Reemplazamos el contenido de la novena columna por "GT"
    }
    print
}' sorted_merged_short.vcf > final_merged_SNP.vcf

#Tras asegurarnos de que el nombre de las muestras en el header de final_merged_SNP.vcf coincide con el nombre de las muestras en fenotipo.csv, lanzamos el análisis GWAS.  

#Activamos nuestro usuario de mamba para poder acceder a vcf2gwas
mamba init
mamba activate ireneggwas
vcf2gwas -v final_merged_SNP.vcf -pf fenotipo.csv -ap -lmm -o output_gwas
#-v se usa para especificar el vcf que se va a analizar. Debe tener los SNP de las muestras fusionados y formateado
#-pf se usa para especificar el archivo de asociación de las muestras con sus fenotipos
#-ap indica que se deben usar todos los fenotipos apra hacer la asociación
#-lmm indica que se debe hacer una asociación usando un modelo lineal mixto 
#-o se usa para indicar el nombre del directorio en el que se va a alamacenar el output
