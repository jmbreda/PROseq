#!/bin/bash


Time=(CT00 CT04 CT08 CT12 CT16 CT20 CT24 CT28 CT32 CT36 CT40 CT44)
N=${#Time[@]}


outfile="resources/sample_to_fastq_Run1.txt"
if [ -f $outfile ]
then
    rm $outfile
fi

for i in $(seq 0 $((N-1)))
do
    printf "%s\t" "${Time[$i]}" >> $outfile

    mapfile -t Files < <( find resources/fastq/PRO_SEQ_CT*_S$((i + 1))_R1_001.fastq.gz )
    n=${#Files[@]}
    for j in $(seq 0 $((n-1)))
    do
        printf "%s" "${Files[$j]}" >> $outfile
        if [ $j -lt $((n-1)) ]
        then
            printf "," >> $outfile
        fi
    done
    printf "\n" >> $outfile
done


outfile="resources/sample_to_fastq_Run2.txt"
if [ -f $outfile ]
then
    rm $outfile
fi

for i in $(seq 0 $((N-1)))
do
    printf "%s\t" "${Time[$i]}" >> $outfile
    for r in $(seq 1 2)
    do
        mapfile -t Files < <( find resources/fastq/18035FL-105-*_S$((i + 1))_*_R"${r}"_001.fastq.gz )
        n=${#Files[@]}
        for j in $(seq 0 $((n-1)))
        do
            printf "%s" "${Files[$j]}" >> $outfile
            if [ $j -lt $((n-1)) ]
            then
                printf "," >> $outfile
            fi
        done

        if [ $r -eq 1 ]
        then
            printf " " >> $outfile
        else
            printf "\n" >> $outfile
        fi
    done
done
