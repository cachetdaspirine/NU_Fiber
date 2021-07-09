#!/bin/bash

for SimNum in {0..9}
do
	sbatch Submit.pbs $SimNum
done
