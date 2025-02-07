#!/bin/bash
#SBATCH --job-name=matlab_compile
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=cidbn
#SBATCH -A cidbn_legacy

# Usage check
if [ $# -lt 2 ]; then
    echo "Usage: sbatch $0 <script_folder> <script_name.m> [MATLAB_version] [log_file]"
    exit 1
fi

# Define input MATLAB script and its directory
SCRIPT_FOLDER=$1
SCRIPT_NAME_WITH_EXT=$2
SCRIPT_NAME=$(basename "$SCRIPT_NAME_WITH_EXT" .m)

# Optional MATLAB version
MATLAB_VERSION=${3:-"default"}  # Default to system MATLAB version if none is provided

# Optional log file (default: "compile_output")
LOG_FILE=${4:-"compile_output"}
#SBATCH --output=${LOG_FILE}

# Check if script folder exists
if [ ! -d "$SCRIPT_FOLDER" ]; then
    echo "Error: Specified folder '$SCRIPT_FOLDER' does not exist."
    exit 1
fi

# Check if MATLAB script exists
SCRIPT_PATH="$SCRIPT_FOLDER/$SCRIPT_NAME_WITH_EXT"
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: MATLAB script '$SCRIPT_NAME_WITH_EXT' not found in '$SCRIPT_FOLDER'."
    exit 1
fi

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

# Compile MATLAB script in the specified folder
echo "Compiling $SCRIPT_NAME_WITH_EXT with MATLAB version $MATLAB_VERSION..." | tee -a "${LOG_FILE}"
cd "$SCRIPT_FOLDER" || exit 1
mcc -m -R -nodisplay -R -singleCompThread -o "$SCRIPT_NAME" "$SCRIPT_NAME_WITH_EXT"

# Print message after successful compilation and write it to the log file
echo "Compiling $SCRIPT_NAME finished!" | tee -a "${LOG_FILE}"
echo "Compilation complete. Output stored in $SCRIPT_FOLDER/" | tee -a "${LOG_FILE}"

