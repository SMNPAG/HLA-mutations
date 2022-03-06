#!/bin/bash
#@seb2022 sebastien.hergalant@inserm.fr

## arg1: dir where sample(s) is(are) - first dot (.) in file name identifies the sample label - this label should be shared between fastq[.gz] and unified report files [*.report-unified.txt]
## arg2: nb threads = how many samples at a time?
## arg3: WGS or TARGETED experiment?

# first check the args
if [[ "$#" -lt 3 ]]; then
    echo "usage: $(basename $0) sample_dir [abs|rel path] nb_parallel_jobs [0 < int < 30] experiment_type [WGS|TARGETED]"
    echo "sample_dir must exist and contain patient .report-unified.txt files along with corresponding fastq[.gz] samples, in paired-end mode"
	echo "nb_parallel_jobs must be a valid integer - within range - and will be used to fork the analysis this number of times"
	echo "experiment_type must be a string among WGS or TARGETED (case insensitive) and represent the origin of the the samples found in sample_dir. Are they coming from WGS (Whole Genome Sequecing) or TARGETED (targeted exome) experiments?"
    exit 1
fi
if [[ ! -d "$1" ]] || [ ! -w "$1" ] ; then echo "$(basename $0): $1 (first argument) is not a suitable directory for fetching samples and writing results" ; exit 1 ; fi
# check if 1 <= $2 <= 29
intregexp='^([1-9]|([1-2][0-9]))$'
if [[ $2 =~ $intregexp ]]; then nb_threads=$2 ; else nb_threads=1 ; echo "$(basename $0): $2 (second argument) is not a valid number. Will assume nb_parallel_jobs = 1" ; fi
# then check if file HLAmutAA_pipeline.conf exists and is sourceable without incident
if [[ ! -r "${0%.*}.conf" ]]; then echo "$(basename $0): good call but config file ${0%.*}.conf must exist and be readable" ; exit 1 ; fi
source "${0%.*}.conf"
if [[ $? -ne 0 ]]; then echo "$(basename $0): there went something wrong when sourcing the config file. Please fix" ; exit 1 ; fi
# check if $3 == TARGETED|WGS - cannot check this before having read the conf file
if [[ ${3^^} != "WGS" ]] && [[ ${3^^} != "TARGETED" ]]; then  echo "$(basename $0): $3 (third argument) is not a suitable experiment type. Valid names are WGS or TARGETED (case insensitive)" ; exit 1 ; fi
export experiment=${3^^}
if [[ $experiment == "TARGETED" ]]; then export VSparamsA=${TARGETEDparamsA} ; export VSparamsF=${TARGETEDparamsF} ; fi
if [[ $experiment == "WGS" ]]; then export VSparamsA=${WGSparamsA} ; export VSparamsF=${WGSparamsF} ; fi

# now all is well. If there is still some problems, the pipeline will throw errors and exit, on a sample by sample basis
export sampleDir=$(echo $1 |sed 's/\/$//') # remove trailing / otherwise some steps will not work properly

## the main "pipeline" function that contains every step of the analysis for a given sample
function exec_sample(){
   # exit this subshell when any command of the function fails
   set -e
   # keep track of the last executed command and print something to help debug before exiting
   trap 'last_command=$current_command; current_command=$BASH_COMMAND; echo "\"${last_command}\" command filed with exit code $?."' ERR

   reportUnified=$1
    
   # dir and prefix variables
   thisSampleDir="${sampleDir}/$(basename ${reportUnified%%.*})"
   samplePrefix="${thisSampleDir}/$(basename ${reportUnified%%.*})"
   [[ -d ${thisSampleDir} ]] && echo "Skipping patient/sample $(basename ${reportUnified%%.*}) as directory ${thisSampleDir} already exists" && continue
   
   (
   # logs for each sample begin here
   echo -ne "****************************************\n\t${experiment} SAMPLE $(basename ${reportUnified%%.*})\n****************************************\n" |tee >(cat >&2)
   echo "*** $(basename ${reportUnified%%.*}) *** Creating directory ${thisSampleDir} for handling this sample $(basename ${reportUnified%%.*}) outputs" |tee >(cat >&2)
   mkdir -p ${thisSampleDir}
   
   ## fasta-parser-OHCC :: to solve the main classes call problem, fix the manifest file
   echo "*** $(basename ${reportUnified%%.*}) *** Extracting fasta sequences for this sample" |tee >(cat >&2)
   java -cp ${fastaParser} us.oh.cc.genomic.FastaParser lociPatientFile=${reportUnified} fastaFile=${hlaDir}/hla_gen.fasta fastaAdditionalFile=${hlaDir}/hla_nuc.fasta >&2
   # no way to control output filenames, so rename by default
   mv ${thisSampleDir}_filtered_allele_0.fasta ${samplePrefix}.allele0.fasta
   mv ${thisSampleDir}_summary\ 0.txt ${samplePrefix}.summary0.txt

   ## novocraft :: index + alignment
   echo "*** $(basename ${reportUnified%%.*}) *** Indexing and aligning with Novocraft" |tee >(cat >&2)
   ${novocraft}/novoindex ${samplePrefix}.novoindex ${samplePrefix}.allele0.fasta >&2
   ${novocraft}/novoalign -t 70 -R 0 -g20 -x 3 -r All 999 -d ${samplePrefix}.novoindex -f ${thisSampleDir}.*R1*.fastq.gz ${thisSampleDir}.*R2*.fastq.gz > ${samplePrefix}.aln.bam ## something to be done here with 1/2 R1/R2 and for fastq or fastq.gz ...

   echo "*** $(basename ${reportUnified%%.*}) *** Sorting, marking duplicates, annotating and computing read depth in the alignment file (bam)" |tee >(cat >&2)
   ## samtools :: prefix sort temp files for parallel use
   ${samtools} sort -m 1G -@ 3 -T ${reportUnified%%.*}.sort.temp -O bam -o ${samplePrefix}.sort.bam ${samplePrefix}.aln.bam
   ## picard :: mark duplicates then annotate bam - no metrics M=/dev/null
   java -Xmx5g -jar ${picard} MarkDuplicates ASSUME_SORTED=true I=${samplePrefix}.sort.bam O=${samplePrefix}.dupl.bam M=/dev/null
   java -jar ${picard} AddOrReplaceReadGroups I=${samplePrefix}.dupl.bam O=${samplePrefix}.group.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$(basename ${reportUnified%%.*})
   ## samtools :: compute alignment depth
   ${samtools} depth ${samplePrefix}.group.bam |awk '{sum+=$3; sumsq+=$3*$3} END { print "Average =",sum/NR; print "Stdev =",sqrt(sumsq/NR - (sum/NR)**2)}' > ${samplePrefix}.group.depth

   ## varscan :: germline mode - SNPs and INDELS, then filter each one of them
   echo "*** $(basename ${reportUnified%%.*}) *** Run VarScan in germline mode for SNPs and INDELs separately (with pileup from samtools0.1.16)" |tee >(cat >&2)
   ## samtools-0.1.16 old pileup version - creates .fai index by default
   ${samtools0116} pileup -f ${samplePrefix}.allele0.fasta ${samplePrefix}.group.bam > ${samplePrefix}.group.pileup
   java -jar ${varscan} pileup2indel ${samplePrefix}.group.pileup ${VSparamsA} > ${samplePrefix}.INDELs
   java -jar ${varscan} pileup2snp ${samplePrefix}.group.pileup ${VSparamsA} > ${samplePrefix}.SNPs
   java -jar ${varscan} filter ${samplePrefix}.SNPs ${VSparamsF} --indel-file ${samplePrefix}.INDELs --output-file ${samplePrefix}.SNPs.filtered
   java -jar ${varscan} filter ${samplePrefix}.INDELs ${VSparamsF} --output-file ${samplePrefix}.INDELs.filtered

   ## fasta-parser :: filter files according to HLA db
   echo "*** $(basename ${reportUnified%%.*}) *** Filter SNPs and INDELs according to the HLA reference database and add topographical locations" |tee >(cat >&2)
   java -cp ${fastaParser} us.oh.cc.genomic.mutations.EnrichVariantCallerFile mutationFile=${samplePrefix}.SNPs.filtered hlaDatFile=${hlaDir}/hla.dat >&2
   mv ${samplePrefix}.txt ${samplePrefix}.SNPs.regions
   java -cp ${fastaParser} us.oh.cc.genomic.mutations.EnrichVariantCallerFile mutationFile=${samplePrefix}.INDELs.filtered hlaDatFile=${hlaDir}/hla.dat >&2
   mv ${samplePrefix}.txt ${samplePrefix}.INDELs.regions
   java -cp ${fastaParser} us.oh.cc.genomic.mutations.MutationAlignmentFilter mutationFile=${samplePrefix}.SNPs.regions alignementFileDirPath=${hlaDir}/alignments/ hlaDatFile=${hlaDir}/hla.dat numberDigitToCheck=0 >&2
   mv ${samplePrefix}_filtered.txt ${samplePrefix}.SNPs.mutAlignFilter
   java -cp ${fastaParser} us.oh.cc.genomic.mutations.MutationInDelFilter mutationFile=${samplePrefix}.INDELs.regions alignementFileDirPath=${hlaDir}/alignments/ hlaDatFile=${hlaDir}/hla.dat numberDigitToCheck=0 >&2
   mv ${samplePrefix}_filtered.txt ${samplePrefix}.INDELs.mutInDelFilter
   
   echo "*** $(basename ${reportUnified%%.*}) *** Cleaning up" |tee >(cat >&2)
   # some cleanup
   rm -f ${samplePrefix}.aln.bam ${samplePrefix}.dupl.bam ${samplePrefix}.sort.bam
   
   echo "*** $(basename ${reportUnified%%.*}) *** Finished with this sample" |tee >(cat >&2)
   echo "$(basename ${reportUnified%%.*}) *** Check log file ${thisSampleDir}.log for verbose outputs"
   set +e
   ) 2> ${thisSampleDir}.log # redirects stderr to log file (all parsers, novocraft, samtools, picard and varscan verbose)
}
export -f exec_sample

## let's go - outputs will be in sampleDir, one dir for each sample
echo "$(basename $0): starting execution"
for report in $(ls -1 ${sampleDir}/*.report_unified.txt |sort -R) ; do echo $report ; done |xargs -I{} --max-procs ${nb_threads} bash -c '{
   report={}
   exec_sample "$report"
}'
echo "$(basename $0): nothing more to do"
