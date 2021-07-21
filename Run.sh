#!/bin/bash

for SimNum in {0..99}
do
	sbatch Submit.pbs $SimNum
done
