#!/bin/bash
#SBATCH --job-name=matlab_compile
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --partition=short

# Usage check
if [ $# -lt 1 ]; then
    echo "Usage: sbatch $0 <matlab_script.m> [MATLAB_version] [log_file]"
    exit 1
fi

# Define input MATLAB script and its directory
MATLAB_SCRIPT=$1
SCRIPT_DIR=$(dirname "$(realpath "$MATLAB_SCRIPT")")
SCRIPT_NAME=$(basename "$MATLAB_SCRIPT" .m)

# Optional MATLAB version
MATLAB_VERSION=${2:-"default"}  # Default to system MATLAB version if none is provided

# Optional log file (default: "compile_output")
LOG_FILE=${3:-"compile_output"}
#SBATCH --output=${LOG_FILE}

# Detect available MATLAB versions
AVAILABLE_MODULES=$(module avail matlab 2>&1)

# Load the correct MATLAB module
if [[ "$MATLAB_VERSION" == "default" ]]; then
    module load matlab  # Load system default MATLAB
elif echo "$AVAILABLE_MODULES" | grep -q "matlab/$MATLAB_VERSION"; then
    module load matlab/"$MATLAB_VERSION"
else
    echo "Error: MATLAB version $MATLAB_VERSION not found." >&2
    exit 1
fi

# Compile MATLAB script in the same directory
echo "Compiling $MATLAB_SCRIPT with MATLAB version $MATLAB_VERSION..." | tee -a "${LOG_FILE}"
cd "$SCRIPT_DIR" || exit 1
mcc -m -R -nodisplay -R -singleCompThread -o "$SCRIPT_NAME" "$MATLAB_SCRIPT"

# Print message after successful compilation and write it to the log file
echo "Compiling $SCRIPT_NAME finished!" | tee -a "${LOG_FILE}"

echo "Compilation complete. Output stored in $SCRIPT_DIR/" | tee -a "${LOG_FILE}"
