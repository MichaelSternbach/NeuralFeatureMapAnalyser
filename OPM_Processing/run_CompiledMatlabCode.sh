#!/bin/bash

# Set default values for variables if not provided
MODULE="${MODULE:-matlab-mcr/R2022b}"
WORKDIR="${WORKDIR:-/scratch/users/sternbach1/OPM_DataPipeline/CompiledCode/}"
EXECUTABLE="${EXECUTABLE:-./OPM_DataPipelineHPC_Faster}"
ARGS="${ARGS:-ferret 1 /home/uni08/cidbn1/ /scratch/users/sternbach1/OPM_DataPipeline/ResultDataTestNewAnalysis/ 1 10 100}"

# Load the specified module
module load "$MODULE"

# Change to the specified working directory
cd "$WORKDIR"

# Run the executable with the provided arguments
"$EXECUTABLE" $ARGS
