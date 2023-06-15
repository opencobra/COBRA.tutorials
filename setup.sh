#!/bin/bash

DEST_REPO=$1
DEST_REPO_TOKEN=$2
FILE_PATH=$3

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

# Extract the directory of the first file
TARGET_DIR="stable/tutorials/$(dirname "$FILE_PATH")"

# Create the target directory in the destination repository if it doesn't exist
mkdir -p $TARGET_DIR

# Delete all contents in the target directory
rm -rf $TARGET_DIR/*


# Add, commit and push the files to the destination repository
git add .
git commit -m "Sync files from source repo" || echo "No changes to commit"
git push origin gh-pages
