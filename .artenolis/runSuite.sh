#!/bin/bash

# save the cloned directory
tutorialsClonedPath=$(pwd)/../
cd $tutorialsClonedPath
mv linux tutorials
mkdir linux && cd linux

# clone the cobratoolbox
git clone --depth=1 --no-single-branch https://github.com/opencobra/cobratoolbox.git cobratoolbox
cd cobratoolbox

# checkout the branch on cobratoolbox (default: develop)
git checkout ci-tutorials-new

# Remove the submodule entry from .git/config
git submodule deinit -f tutorials

# Remove the submodule directory from the superproject's .git/modules directory
rm -rf .git/modules/tutorials

# Remove the entry in .gitmodules and remove the submodule directory located at path/to/submodule
git rm -f tutorials

# commit the removed submodule
git commit -m "Removed tutorials submodule"

# initialize the submodules
git submodule update --init --depth=1 --remote --no-fetch

# move the cloned tutorials folder to the cobratoolbox directory
cd $tutorialsClonedPath
mv tutorials linux/cobratoolbox/.
cd linux/cobratoolbox

# launch the script
bash .artenolis/runTutorials.sh