## GET FASTA FILE FOR GENCODE ANNOTATIONS
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
gunzip gencode.v38.transcripts.fa.gz

salmon index --gencode -t gencode.v38.transcripts.fa -i human_hg38_index -k 21

numBootstraps=50
numThreads=8

input="SPREADSHEET_WITH_FASTQ_PATHS.txt"
cat ${input} | while read line
    do
    echo "Processing file ${line}"
    name=$(basename $f .fastq.gz)
    ### MAKE A NAME FOR OUTNAME
    
    salmon quant \
      -i human_hg38_index \
      -1 ${read1} \
      -2 ${read2} \
      --validateMappings \
      --useEM \
      --numBootstraps ${numBootstraps} \
      --seqBias \
      --gcBias \
      -p ${numThreads} \
      -o ${outName}
    
done 