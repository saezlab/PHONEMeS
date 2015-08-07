#!/bin/bash

j=1

genM=p15gen$j

for i in {1..50}

	do
	bsub  -J "$genM" -q research-rh6 -M 200 -R "rusage[mem=200]" Rscript scriptGxopt_50models_15.R 1 $i $j
	
done

genRes=p15res$j

bsub  -J "$genRes" -w "done($genM)" -q research-rh6 -M 1500 -R "rusage[mem=1500]" Rscript processGx_15.R 1 $j

for j in {2..50}

do
	genM=p15gen$j
	for i in {1..50}

		do
		bsub  -J "$genM" -w "done($genRes)" -o "p15out" -q research-rh6 -M 200 -R "rusage[mem=200]" Rscript scriptGxopt_50models_15.R 1 $i $j
	
	done
	
	genRes=p15res$j
	
	bsub  -J "$genRes" -w "done($genM)" -q research-rh6 -M 1500 -R "rusage[mem=1500]" Rscript processGx_15.R 1 $j
	
done	

bsub  -J "p15import" -w "done(p15res50)" -q research-rh6 -M 204800 -R "rusage[mem=204800]" Rscript import_15.R