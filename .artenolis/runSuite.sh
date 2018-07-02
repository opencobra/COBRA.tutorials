#!/bin/bash

tutorialsClonedPath=$(pwd)/../
cd $tutorialsClonedPath
mv linux tutorials
mkdir linux
cd linux
git clone --depth=1 --no-single-branch https://github.com/opencobra/cobratoolbox.git cobratoolbox
cd cobratoolbox
git checkout ci-tutorials-new # checkout the branch on cobratoolbox (default: develop)

# Remove the submodule entry from .git/config
git submodule deinit -f tutorials

# Remove the submodule directory from the superproject's .git/modules directory
rm -rf .git/modules/tutorials

# Remove the entry in .gitmodules and remove the submodule directory located at path/to/submodule
git rm -f tutorials

git commit -m "Removed tutorials submodule"

git submodule update --init --depth=1 --remote --no-fetch


cd $tutorialsClonedPath # move out of the cobratoolbox directory
mv tutorials linux/cobratoolbox/.
cd linux/cobratoolbox
bash .artenolis/runTutorials.sh