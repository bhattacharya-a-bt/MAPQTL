numBootstraps=50
numThreads=8

input="salmon_rnaseq_filenames_ELGAN1.txt"
cat ${input} | while read line
    do
    echo "Processing file ${line}"
    basename=`echo $line | sed 's/^.*\(LS.*-RNA\).*$/\1/'`
    outName=${basename}_quant
    

sbatch -p general -N 1 --mem=96g -n 1 -t 05-12:00:00 --wrap="salmon quant \
      -i human_hg38_index \
      -l A \
      -r ${line} \
      --validateMappings \
      --useEM \
      --numBootstraps ${numBootstraps} \
      --seqBias \
      --gcBias \
      -p ${numThreads} \
      -o ${outName}"
    
done 