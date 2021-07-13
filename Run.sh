#!/bin/bash

for SimNum in {0..19}
do
	sbatch Submit.pbs $SimNum
done
