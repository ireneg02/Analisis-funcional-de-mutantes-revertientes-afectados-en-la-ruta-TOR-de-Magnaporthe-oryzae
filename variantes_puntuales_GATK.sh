#!/usr/bin/env bash
#SBATCH --mem=200GB
#SBATCH --job-name variantes_puntuales
#SBATCH --partition=medium
#SBATCH --time=2-00:00:00          # total run time limit in HH:MM:SS

#Cargamos los módulos y el programa:
module use /beegfs/easybuild/CentOS/7.6.1810/Skylake/modules/all
module use /beegfs/easybuild/common/modules/all
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
module load picard/2.26.10-Java-15
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load BCFtools/1.9-GCC-8.2.0-2.31.1

#Antes de detectar las variantes puntuales, obtenemos el archivo BAM a partir del SAM, lo ordenamos y lo indexamos.
#En los mapeados de BWA faltaban los grupos de lecturas, por lo que los añadimos manualmente antes de procesar los archivos para la detección de variantes
cd /home/ireneg/DATOS_DNA-SEQ/procesamiento/ #Nos desplazamos a la carpeta con las carpetas de las muestras
for folder in $(ls -d */); do #Recorremos las carpetas de las muestras
        cd $folder #Entramos en la carpeta
        file="${folder:0:18}" #Obtenemos el nombre de la muestra
        
        #Primero procesamos los archivos de BWA
        samtools view -Sb -o BWA/bwa_map_${file}_merge_PE_SE_Guy11.bam BWA/bwa_map_${file}_merge_PE_SE_Guy11.sam #Convertimos el archivo SAM a BAM
        samtools sort -O bam -o BWA/sorted_bwa_map_${file}_merge_PE_SE_Guy11.bam BWA/bwa_map_${file}_merge_PE_SE_Guy11.bam #Ordenamos el archivo BAM
        java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=BWA/sorted_bwa_map_${file}_merge_PE_SE_Guy11.bam O=BWA/sorted_bwa_map_${file}_merge_PE_SE_Guy11_group.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20 #Añadimos los grupos a las lecturas
        samtools index BWA/sorted_bwa_map_${file}_merge_PE_SE_Guy11_group.bam #Indexamos el archivo resultante
        
        #Identificamos variantes
        gatk HaplotypeCaller --input "BWA/sorted_bwa_map_${file}_merge_PE_SE_group.bam" --output "GATK_BWA/GATK_BWA_${file}.vcf" --reference "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" #Identificamos las variantes de los archivos de mapeado BWA respecto al genoma de referencia Guy11
        #Filtramos las variantes detectadas
        gatk SelectVariants -V "GATK_BWA/GATK_BWA_${file}.vcf" --select-type-to-include SNP -O "GATK_BWA/SNP_GATK_BWA_${file}.vcf" #Obtenemos un archivo con únicamente los SNP detectados
        gatk VariantFiltration -R "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" -V GATK_BWA/SNP_GATK_BWA_${file}.vcf -O GATK_BWA/SNP_qualfilter_all_BWA_${file}.vcf --filter-name One --filter-expression "QD < 2.0" --filter-name Two --filter-expression "AF < 1.0" --filter-name Three --filter-expression "FS > 60.0" --filter-name Four --filter-expression "MQ < 40.0" --filter-name Five --filter-expression  "MQRankSum < -12.5" --filter-name Six --filter-expression "ReadPosRankSum < -8.0" #Aplicamos filtros para eliminar los SNP de baja calidad
        bcftools view -i 'FILTER="PASS"' GATK_BWA/SNP_qualfilter_all_BWA_${file}.vcf -o GATK_BWA/SNP_qualfilter_BWA_${file}.vcf #Obtenemos un archivo con los SNP que han pasado los filtros
        cp GATK_BWA/SNP_qualfilter_BWA_${file}.vcf ../gwas #Copiamos el archivo resultante en una carpeta en la que posteriormente realizaremos el análisis GWAS
        
        gatk SelectVariants -V "GATK_BWA/GATK_BWA_${file}.vcf" --select-type-to-include INDEL -O "GATK_BWA/INDEL_GATK_BWA_${file}.vcf" #Obtenemos un archivo con únicamente los indel detectados
        gatk VariantFiltration -R "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" -V GATK_BWA/INDEL_GATK_BWA_${file}.vcf -O GATK_BWA/INDEL_qualfilter_all_BWA_${file}.vcf --filter-name One --filter-expression "QD < 2.0" --filter-name Two --filter-expression "AF < 1.0" --filter-name Three --filter-expression "FS > 200.0" --filter-name Four --filter-expression "ReadPosRankSum < -20.0" #Aplicamos filtros para eliminar los indel de baja calidad
        bcftools view -i 'FILTER="PASS"' GATK_BWA/INDEL_qualfilter_all_BWA_${file}.vcf -o GATK_BWA/INDEL_qualfilter_BWA_${file}.vcf #Obtenemos un archivo con los indel que han pasado los filtros
        
        
        #A continuación, procesamos los archivos de Segemehl
        samtools view -Sb -o SEGEMEHL/segemehl_map_${file}_merge_PE_SE_Guy11.bam SEGEMEHL/segemehl_map_${file}_merge_PE_SE_Guy11.sam #Convertimos el archivo SAM a BAM
        samtools sort -O bam -o SEGEMEHL/sorted_segemehl_map_${file}_merge_PE_SE_Guy11.bam SEGEMEHL/segemehl_map_${file}_merge_PE_SE_Guy11.bam #Ordenamos el archivo BAM 
        #Identificamos las variantes
        gatk HaplotypeCaller --input "SEGEMEHL/sorted_segemehl_map_${file}_merge_PE_SE.bam" --output "GATK_SEGEMEHL/GATK_SEGEMEHL_${file}.vcf" --reference "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" #Identificamos las variantes de los archivos de mapeado Segemehl respecto al genoma de referencia Guy11
        #Filtramos las variantes detectadas
        gatk SelectVariants -V "GATK_SEGEMEHL/GATK_SEGEMEHL_${file}.vcf" --select-type-to-include SNP -O "GATK_SEGEMEHL/SNP_GATK_SEGEMEHL_${file}.vcf" #Obtenemos un archivo con únicamente los SNP detectados
        gatk VariantFiltration -R "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" -V GATK_SEGEMEHL/SNP_GATK_SEGEMEHL_${file}.vcf -O GATK_SEGEMEHL/SNP_qualfilter_all_SEGEMEHL_${file}.vcf --filter-name One --filter-expression "QD < 2.0" --filter-name Two --filter-expression "AF < 1.0" --filter-name Three --filter-expression "FS > 60.0" --filter-name Four --filter-expression "MQ < 40.0" --filter-name Five --filter-expression  "MQRankSum < -12.5" --filter-name Six --filter-expression "ReadPosRankSum < -8.0" #Aplicamos filtros para eliminar los SNP de baja calidad
        bcftools view -i 'FILTER="PASS"' GATK_SEGEMEHL/SNP_qualfilter_all_SEGEMEHL_${file}.vcf -o GATK_SEGEMEHL/SNP_qualfilter_SEGEMEHL_${file}.vcf #Obtenemos un archivo con los SNP que han pasado los filtros
        gatk SelectVariants -V "GATK_SEGEMEHL/GATK_SEGEMEHL_${file}.vcf" --select-type-to-include INDEL -O "GATK_SEGEMEHL/INDEL_GATK_SEGEMEHL_${file}.vcf" #Obtenemos un archivo con únicamente los indel detectados
        gatk VariantFiltration -R "/home/ireneg/DATOS_DNA-SEQ/originales/MagorGY11_1/Mycocosm/Assembly/Genome Assembly (unmasked)/MagorGY11_1_AssemblyScaffolds.fasta" -V GATK_SEGEMEHL/INDEL_GATK_SEGEMEHL_${file}.vcf -O GATK_SEGEMEHL/INDEL_qualfilter_all_SEGEMEHL_${file}.vcf --filter-name One --filter-expression "QD < 2.0" --filter-name Two --filter-expression "AF < 1.0" --filter-name Three --filter-expression "FS > 200.0" --filter-name Four --filter-expression "ReadPosRankSum < -20.0" #Aplicamos filtros para eliminar los indel de baja calidad
        bcftools view -i 'FILTER="PASS"' GATK_SEGEMEHL/INDEL_qualfilter_all_SEGEMEHL_${file}.vcf -o GATK_SEGEMEHL/INDEL_qualfilter_SEGEMEHL_${file}.vcf #Obtenemos un archivo con los indel que han pasado los filtros
        
        cd .. #Salimos de la carpeta de la muestra
done



#En el comando de conversión de SAM a BAM usando samtools view: 
#-Sb Convierte de formato SAM a BAM.
#-o permite especificar a continuación el nombre del archivo de salida en formato BAM.
#Finalmente se aporta el archivo SAM a convertir

#En el comando para ordenar el archivo BAM con samtools sort
#-O bam: Especifica el formato de salida como BAM.
#-o permite especificar a continuación el nombre del archivo de salida
#Finalmente se aporta el archivo BAM a ordenar

#En el comando para añadir grupos de lectura con Picard:
#I permite indicar el archivo de entrada
#O permite indicar el nombre del archivo de salida.
#RGID=4, RGLB=lib1, RGPL=ILLUMINA, RGPU=unit1, RGSM=20 son los parámetros que definen los grupos de lectura. 

#En el comando para indexar el archivo BAM con samtools index, simplemente se indica el archivo a indexar

#Para la identificación de variantes se usa gatk HaplotypeCaller:
#--input a el archivo BAM de entrada
#--output especifica el archivo VCF de salida con las variantes
#--reference indica el archivo FASTA de referencia

#gatk SelectVariants permite hacer una selección de las variantes
#-V indica el archivo VCF de entrada
#--select-type-to-include permite indicar a continuación el tipo de variante que queremos seleccionar
#-O especifica el nombre del archivo VCF que se va a generar con las variantes seleccionadas

#gatk VariantFiltration permite hacer un filtrado de variantes en función a la calidad
#-R  indica el archivo FASTA de referencia
#-V indica el archivo VCF de entrada
#-O especifica el nombre del archivo VCF que se va a generar, en el que se indica si las varainates pasan o no los filtros
#--filter-name y --filter-expression son parámetros para definir los criterios de filtrado (QD, AF, FS, MQ, MQRankSum, ReadPosRankSum).

#Para seleccionar las variantes que pasaron los filtros usamos bcftools view 
#-i 'FILTER="PASS"' indica que se seleccionan solo las variantes que pasaron todos los filtros
#Luego se indica el archivo de entrada
#-o indica el nombre del archivo que se genera con las variantes que pasaron los filtros. 
