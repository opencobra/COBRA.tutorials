#!/bin/bash

DEST_REPO=$1
FILE_PATH=$2

echo "Starting revert script execution..."

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

# Reset the branch to the previous commit
echo "Resetting to the previous commit..."
git reset --hard HEAD~1

# Force push the reset branch
echo "Force pushing the revert..."
git push origin gh-pages --force

echo "Revert script execution completed."
