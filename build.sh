#!/bin/bash

DEST_REPO=$1
FILE_PATH=$2

echo "Starting script execution..."

# Convert mlx to html and pdf using MATLAB command
# Note: Ensure that MATLAB is installed and added to the PATH
echo "Converting .mlx file to .html and .pdf..."
ABSOLUTE_FILE_PATH=$(realpath "$FILE_PATH")
HTML_FILE_PATH=$(echo "$ABSOLUTE_FILE_PATH" | sed 's/.mlx/.html/g')
PDF_FILE_PATH=$(echo "$ABSOLUTE_FILE_PATH" | sed 's/.mlx/.pdf/g')
echo "Absolute file path: $ABSOLUTE_FILE_PATH"
echo "HTML file path: $HTML_FILE_PATH"
echo "PDF file path: $PDF_FILE_PATH"

# Setup virtual frame buffer
export DISPLAY=:100
echo "Starting virtual frame buffer..."
Xvfb -ac :100 -screen 0 1280x1024x24 > /dev/null &

# Run MATLAB command to convert .mlx to .html and .pdf
echo "Running MATLAB conversion command..."
/usr/local/MATLAB/R2024a/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$ABSOLUTE_FILE_PATH', '$HTML_FILE_PATH')"
echo "Converting to PDF..."
/usr/local/MATLAB/R2024a/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$ABSOLUTE_FILE_PATH', '$PDF_FILE_PATH')"

# Clone the destination repository
echo "Cloning the destination repository: https://github.com/$DEST_REPO.git"
git clone --depth 1 git@github.com:$DEST_REPO.git

# Split the destination repository into owner and name
IFS='/' read -ra ADDR <<< "$DEST_REPO"
DEST_REPO_OWNER=${ADDR[0]}
DEST_REPO_NAME=${ADDR[1]}
echo "Destination repository owner: $DEST_REPO_OWNER"
echo "Destination repository name: $DEST_REPO_NAME"

# Change to the destination repository directory
echo "Changing to the destination repository directory: $DEST_REPO_NAME"
cd $DEST_REPO_NAME

# Set up git config
echo "Setting up git config..."
git config user.name "github-actions[bot]"
git config user.email "github-actions[bot]@users.noreply.github.com"

# Create the target directory in the destination repository
TARGET_DIR="docs/tutorials/$(dirname "$FILE_PATH")"
echo "Creating the target directory: $TARGET_DIR"
mkdir -p "$TARGET_DIR"

# Copy the HTML, PDF, and .mlx files to the target directory in the destination repository
echo "Copying the HTML, PDF, and .mlx files to the target directory..."
cp "$HTML_FILE_PATH" "$TARGET_DIR/"
cp "$PDF_FILE_PATH" "$TARGET_DIR/"
cp "$ABSOLUTE_FILE_PATH" "$TARGET_DIR/"

# Add, commit, and push the changes
echo "Adding changes to git..."
git add .
echo "Committing changes..."
git commit -m "Sync files from source repo" || echo "No changes to commit"
echo "Pushing changes to docs folder in master branch..."
git push origin master

echo "Script execution completed."
