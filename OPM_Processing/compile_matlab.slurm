#!/bin/bash
#SBATCH --job-name=matlab_compile
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=cidbn
#SBATCH -A cidbn_legacy

# Usage check
if [ $# -lt 2 ]; then
    echo "Usage: sbatch $0 <script_folder> <script_name.m> [MATLAB_version]"
    exit 1
fi

# Inputs
SCRIPT_FOLDER=$1
SCRIPT_NAME=$2
MATLAB_VERSION=$3

# Change to the specified folder
cd "$SCRIPT_FOLDER" || exit 1

# Load MATLAB
module load matlab/${MATLAB_VERSION} 2>/dev/null || module load matlab

# Compile the script
mcc -m -R -nodisplay -R -singleCompThread "$SCRIPT_NAME"

# Echo completion
echo "Compiling $SCRIPT_NAME finished!"

