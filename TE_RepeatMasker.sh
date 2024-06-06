#Esta herramienta se ha utilizado desde Galaxy, pero la transcripción del comando es:

# Definir variables
INPUT_FASTA="/path/to/input_fasta_file.fasta"
OUTPUT_DIR="/path/to/output_directory"
MASKED_OUTPUT="${OUTPUT_DIR}/masked_output.dat"
OUT_FILE="${OUTPUT_DIR}/output_file.dat"
TBL_FILE="${OUTPUT_DIR}/table_file.dat"
GFF_FILE="${OUTPUT_DIR}/output_file.gff"
ALIGN_FILE="${OUTPUT_DIR}/align_file.dat"
CAT_OUTPUT="${OUTPUT_DIR}/cat_output.dat"

# Definir el path de RepeatMasker
RM_PATH=$(which RepeatMasker)
if [ -z "$RM_PATH" ]; then
  echo "Failed to find RepeatMasker in PATH ($PATH)" >&2
  exit 1
fi

# Definir el path de la librería que se va a usar
if [ -z "$RM_LIB_PATH" ]; then
  RM_LIB_PATH=$(dirname "$RM_PATH")/../share/RepeatMasker/Libraries
fi

# Crear un link simbólico al archivo fasta de entrada 
ln -s "${INPUT_FASTA}" rm_input.fasta

# Correr RepeatMasker
RepeatMasker -dir "$(pwd)" -libdir "$RM_LIB_PATH" -lib "${INPUT_FASTA}" \
  -cutoff '225' -parallel "${GALAXY_SLOTS:-1}" -gff -excln -nolow -s \
  -frag 40000 -ali rm_input.fasta

# Mover and procesar los archivos de output
mv rm_input.fasta.masked "${MASKED_OUTPUT}"
sed -E 's/^ *//; s/ *$//; s/\+ //; s/ +/\t/g; 1,2c SW score\t% div.\t% del.\t% ins.\tquery sequence\tpos in query: begin\tend\t(left)\trepeat\tclass/family\tpos in repeat: begin\tend\t(left)\tID' \
  rm_input.fasta.out > "${OUT_FILE}"
mv rm_input.fasta.tbl "${TBL_FILE}"
mv rm_input.fasta.out.gff "${GFF_FILE}"
mv rm_input.fasta.align "${ALIGN_FILE}"

# Manipular el archivo .cat
if [ -f 'rm_input.fasta.cat.gz' ]; then
  zcat 'rm_input.fasta.cat.gz' > "${CAT_OUTPUT}"
else
  mv rm_input.fasta.cat "${CAT_OUTPUT}"
fi
