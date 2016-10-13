#!/bin/bash

j=1

genM=p${1}gen$j

for i in {1..50}

	do
	bsub -P jrc_combine -J "$genM" -M 200 -R "rusage[mem=200]" Rscript scriptGxopt_50models.R $1 $i $j
	
done

genRes=p${1}res$j
bsub -P jrc_combine -J "$genRes" -w "done($genM)" -M 1500 -R "rusage[mem=1500]" Rscript processGx.R $1 $j

for j in {2..50}

do
	genM=p${1}gen$j
	for i in {1..50}

		do
		bsub -P jrc_combine -J "$genM" -w "done($genRes)" -o "p${1}out" -M 200 -R "rusage[mem=200]" Rscript scriptGxopt_50models.R $1 $i $j
	
	done
	
	genRes=p${1}res$j
	
	bsub -P jrc_combine -J "$genRes" -w "done($genM)" -M 1500 -R "rusage[mem=1500]" Rscript processGx.R $1 $j
	
done	

bsub -P jrc_combine -J "p${1}import" -w "done(p${1}res50)" -M 20000 -W 100:00 Rscript import_x.R $1
