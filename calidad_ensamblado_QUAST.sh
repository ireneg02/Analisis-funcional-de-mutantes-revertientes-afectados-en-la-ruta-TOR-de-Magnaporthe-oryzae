#Este análisis se realizó en Galaxy, pero la transcripción de las opciones usadas es:

# Definir variables
INPUT_FASTA="/path/to/input_fasta_file.fasta"
REFERENCE="/path/to/reference_genome.fasta"
OUTPUT_DIR="/path/to/output_directory"
THREADS=${GALAXY_SLOTS:-1}
FINAL_OUTPUT_DIR="/path/to/final_output_directory"

# Correr QUAST
quast --labels 'genome_scf_fasta' \
      -o "${OUTPUT_DIR}" \
      -r "${REFERENCE}" \
      --fungus \
      --min-identity 95.0 \
      --min-contig 500 \
      --split-scaffolds \
      --min-alignment 65 \
      --ambiguity-usage 'one' \
      --ambiguity-score 0.99 \
      --local-mis-size 200 \
      --contig-thresholds '0,1000' \
      --extensive-mis-size 1000 \
      --scaffold-gap-max-size 1000 \
      --unaligned-part-size 500 \
      --x-for-Nx 90 \
      "${INPUT_FASTA}" \
      --threads "${THREADS}"

# Copiar los resultados al directorio final
mkdir -p "${FINAL_OUTPUT_DIR}"
cp "${OUTPUT_DIR}"/*.html "${FINAL_OUTPUT_DIR}"
cp -R "${OUTPUT_DIR}/icarus_viewers" "${FINAL_OUTPUT_DIR}"
