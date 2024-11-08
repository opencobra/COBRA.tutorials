#!/bin/bash

DEST_REPO=$1
FILE_PATH=$2

echo "Starting script execution..."

# Clone the destination repository using SSH
echo "Cloning the destination repository: git@github.com:$DEST_REPO.git"
git clone git@github.com:$DEST_REPO.git

# Split the destination repository into owner and name
IFS='/' read -ra ADDR <<< "$DEST_REPO"
DEST_REPO_OWNER=${ADDR[0]}
DEST_REPO_NAME=${ADDR[1]}
echo "Destination repository owner: $DEST_REPO_OWNER"
echo "Destination repository name: $DEST_REPO_NAME"

# Change the current directory to the destination repository
echo "Changing to the destination repository directory: $DEST_REPO_NAME"
cd $DEST_REPO_NAME

# Set up git config
echo "Setting up git config..."
git config user.name "GitHub Action"
git config user.email "action@github.com"

echo "Checking out gh-pages branch..."
git checkout gh-pages

# Extract the directory of the first file
TARGET_DIR="stable/tutorials/$(dirname "$FILE_PATH")"
echo "Target directory: $TARGET_DIR"

# Check if the target directory is exactly "stable/tutorials"
if [[ "$TARGET_DIR" == "stable/tutorials" ]]; then
    echo "Target directory is 'stable/tutorials'. Do not change."
else
    echo "Target directory is not 'stable/tutorials'. Proceeding with creation."

    # Create the target directory in the destination repository if it doesn't exist
    echo "Creating the target directory if it doesn't exist..."
    mkdir -p $TARGET_DIR
fi

# Add, commit, and push the files to the destination repository
echo "Adding changes to git..."
git add .
echo "Committing changes..."
git commit -m "Sync files from source repo" || echo "No changes to commit"
echo "Pushing changes to gh-pages branch..."
git push origin gh-pages

echo "Setup.sh Script execution completed."
