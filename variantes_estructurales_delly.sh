#Este análisis se realizó en Galaxy, pero la transcripción de las opciones usadas es:

# Definir variables
BAM_FILE="/path/to/input_0.bam"
BAM_INDEX_FILE="/path/to/input_0.bam.bai"
GENOME_FILE="/path/to/genome.fa"
OUTFILE="result.bcf"
VCF_FILE="result.vcf"

# Crear links simbólicos
ln -s "${BAM_FILE}" 'input_0.bam'
ln -s "${BAM_INDEX_FILE}" 'input_0.bam.bai'
ln -s "${GENOME_FILE}" 'genome.fa'

# Correr Delly
delly call --svtype ALL --genome genome.fa --outfile "${OUTFILE}" --map-qual 1 --qual-tra 20 --mad-cutoff 9 --minclip 25 --min-clique-size 2 --minrefsep 25 --maxreadsep 40 --geno-qual 5 'input_0.bam'

# Comprobar si result.bcf existe y convertir a formato VCF
if test -f "${OUTFILE}"; then
  bcftools view "${OUTFILE}" > "${VCF_FILE}"
else
  echo 'No results.'
fi
