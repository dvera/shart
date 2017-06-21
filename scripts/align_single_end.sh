#!/bin/bash

#fastq=$1
#index_file=$2
#prefix=$3
#nThreads=$4

while getopts ":i:t:" opt; do
  case $opt in
    t)
      if [[ "$OPTARG"=="-*" ]]; then
				((OPTIND--))
	    	NTHREADS=1
			else
				NTHREADS=$OPTARG
      fi
      ;;
    i)
      INDEXFILE=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    [?])	print >&2 "Usage: $0 [-i] bwaIndex [-t threads] file1 [file2 fileN ... ]"
		  exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

FASTQFILES=$@

if [[ -f ${INDEXFILE}.bwt ]]; then
  INDEXPREFIX=$INDEXFILE
elif [[ -f $INDEXFILE ]] && [[ $INDEXFILE=~ \.tar\.gz$ ]]
  mkdir -p bwaIndex
  INDEXPATH=$(readlink -f $INDEXFILE)
  tar -C bwaIndex -x -f $INDEXPATH
  INDEXPREFIX=$(readlink -f bwaIndex/*.bwt | sed 's/\.bwt$//g')
else
  echo "index not found"
  exit 2
fi

# run bwa
for f in $FASTQFILES; do
  # unzip file if compressed
  if [[ ! -f $f ]]; then
    echo "fastq file not found"
    exit 1
  elif [[ $f =~ \.gz$ ]]
    fo=${fastq/%\.gz/}
    gunzip -c $f > $fo
    f=$fo
  fi
  
  # create output prefix
  OUTPUT="$(basename $f | sed 's/\.fq$//g' | sed 's/\.fastq$//g').bam"
  
  # check to see input name has same as output name
  #   if [[ $OUTPUT == $f ]]; then
  #     echo "fastq file has inappropriate name, must be a fastq"
  #     exit 1
  #   fi
  
  # align fastq file and run samstats
  bwa mem -t $NTHREADS $INDEXPREFIX $f | samtools view -Shb - > $OUTPUT
  samtools stats $OUTPUT > ${OUTPUT}.samstats
done
 
 
