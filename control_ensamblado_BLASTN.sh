#Este programa se ha ejecutado desde Galaxy, pero la equivalencia de comandos es la siguiente:

# Definr variables
QUERY_FILE="/path/to/query_file.dat"
SUBJECT_FILE="/path/to/subject_file.dat"
OUTPUT_FILE="/path/to/output_directory/output_file.dat"
NUM_THREADS="${GALAXY_SLOTS:-8}"

# Correr blastn
blastn -query "${QUERY_FILE}" \
       -subject "${SUBJECT_FILE}" \
       -task 'blastn' \
       -evalue '0.001' \
       -out "${OUTPUT_FILE}" \
       -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' \
       -num_threads "${NUM_THREADS}"
