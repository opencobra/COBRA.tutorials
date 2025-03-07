#!/bin/bash

DEST_REPO=$1
FILE_PATH=$2

echo "Starting script execution..."

# Convert mlx to html, pdf and m using MATLAB command
# Note: Ensure that MATLAB is installed and added to the PATH
echo "Converting .mlx file to .HTML, .pdf and .m..."
ABSOLUTE_FILE_PATH=$(realpath "$FILE_PATH")
HTML_FILE_PATH=$(echo "$ABSOLUTE_FILE_PATH" | sed 's/.mlx/.html/g')
PDF_FILE_PATH=$(echo "$ABSOLUTE_FILE_PATH" | sed 's/.mlx/.pdf/g')
M_FILE_PATH=$(echo "$ABSOLUTE_FILE_PATH" | sed 's/.mlx/.m/g')
echo "Absolute file path: $ABSOLUTE_FILE_PATH"
echo "HTML file path: $HTML_FILE_PATH"
echo "PDF file path: $PDF_FILE_PATH"
echo "PDF file path: $M_FILE_PATH"

# Run MATLAB command to convert  .mlx file to .HTML, .pdf and .m..."
echo "Running MATLAB conversion command..."
/usr/local/MATLAB/R2024a/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$ABSOLUTE_FILE_PATH', '$HTML_FILE_PATH')"
/usr/local/MATLAB/R2024a/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$ABSOLUTE_FILE_PATH', '$PDF_FILE_PATH')"
/usr/local/MATLAB/R2024a/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$ABSOLUTE_FILE_PATH', '$M_FILE_PATH')"

# Split the destination repository into owner and name
IFS='/' read -ra ADDR <<< "$DEST_REPO"
DEST_REPO_OWNER=${ADDR[0]}
DEST_REPO_NAME=${ADDR[1]}
echo "Destination repository owner: $DEST_REPO_OWNER"
echo "Destination repository name: $DEST_REPO_NAME"

# Change to the destination repository directory
echo "Changing to the destination repository directory: $DEST_REPO_NAME"
cd $DEST_REPO_NAME

# checking out gh-pages branch
git checkout gh-pages

# Create the target directory in the destination repository
TARGET_DIR="stable/tutorials/$(dirname "$FILE_PATH")"
echo "Creating the target directory: $TARGET_DIR"
mkdir -p "$TARGET_DIR"

# Copy the HTML, PDF, mlx and .m files to the target directory in the destination repository
echo "Copying the HTML, PDF, mlx and .m files to the target directory..."
cp "$HTML_FILE_PATH" "$TARGET_DIR/"
cd ../
