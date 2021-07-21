#!/bin/bash

for SimNum in {60..79}
do
	sbatch Submit.pbs $SimNum
done
