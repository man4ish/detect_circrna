input_file="in/reads/jakobi2016_sra_list.txt"

while read -r srr_id; do
  echo "Downloading paired-end reads for $srr_id ..."
  fasterq-dump --split-files "$srr_id"
done < "$input_file"
