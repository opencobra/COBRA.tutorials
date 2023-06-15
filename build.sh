#!/bin/bash

DEST_REPO=$1
DEST_REPO_TOKEN=$2
FILE_PATH=$3

# Convert mlx to html using MATLAB command
# Note: Ensure that MATLAB is installed and added to the PATH
HTML_FILE_PATH=$(echo "$FILE_PATH" | sed 's/.mlx/.html/g')
export DISPLAY=:99
Xvfb -ac :99 -screen 0 1280x1024x24 > /dev/null & /home/aaron/Documents/Matlab/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$FILE_PATH', '$HTML_FILE_PATH')"

# Clone the destination repository
git clone https://AaronBrennan1:$DEST_REPO_TOKEN@github.com/$DEST_REPO.git

# Split the destination repository into owner and name
IFS='/' read -ra ADDR <<< "$DEST_REPO"
DEST_REPO_OWNER=${ADDR[0]}
DEST_REPO_NAME=${ADDR[1]}

# Change the current directory to the destination repository
cd $DEST_REPO_NAME

# Set up git config
git config user.name "GitHub Action"
git config user.email "action@github.com"

git checkout gh-pages

# Create the target directory in the destination repository
mkdir -p "stable/tutorials/$(dirname "$HTML_FILE_PATH")"

# Copy the file to the target directory in the destination repository
cp "../$HTML_FILE_PATH" "stable/tutorials/$HTML_FILE_PATH"

# Add, commit and push the files to the destination repository
git add .
git commit -m "Sync files from source repo" || echo "No changes to commit"
git push origin gh-pages
