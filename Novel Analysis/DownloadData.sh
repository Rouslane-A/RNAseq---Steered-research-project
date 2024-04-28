# loading the module sratoolkit
module load sratoolkit/3.0.0-5fetwpi

# Path to the file containing the list of SRR accessions downloaded from the SRA website from the GEO
accession_list="SRR_Acc_List.txt"

# Loop through each accession in the list, 
while IFS= read accession; do
	prefetch $accession  # Download SRA data
	fasterq-dump --progress --split-files $accession  # Convert to FASTQ
    	rm -rf $accession
done < "$accession_list"
