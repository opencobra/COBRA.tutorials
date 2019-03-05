#!/bin/bash

if [ "$JENKINS" == "True" ]; then
    # save the cloned directory
    tutorialsClonedPath=$(pwd)/../
    cd $tutorialsClonedPath

    # rename the current cloned folder as tutorials (depends on the label of the job)
    mv $ARTENOLIS_SLAVE_NAME_01 tutorials
    mkdir $ARTENOLIS_SLAVE_NAME_01 && cd $ARTENOLIS_SLAVE_NAME_01

    # clone the cobratoolbox
    git clone --depth=1 --no-single-branch https://github.com/opencobra/cobratoolbox.git cobratoolbox
    cd cobratoolbox

    # checkout the branch on cobratoolbox (default development branch: develop)
    git checkout develop

    # remove the submodule entry from .git/config
    git submodule deinit -f tutorials

    # remove the submodule directory from the superproject's .git/modules directory
    rm -rf .git/modules/tutorials

    # remove the entry in .gitmodules and remove the submodule directory located at path/to/submodule
    git rm -f tutorials

    # commit the removed submodule
    git commit -m "[temporary commit] remove tutorials submodule"

    # initialize the submodules
    git submodule update --init --remote --no-fetch

    # move the cloned tutorials folder to the cobratoolbox directory
    cd $tutorialsClonedPath
    mv tutorials $ARTENOLIS_SLAVE_NAME_01/cobratoolbox/.
    cd $ARTENOLIS_SLAVE_NAME_01/cobratoolbox
fi

# define the path to the tutorials
COBRATutorialsPath=$(pwd)

# gather the list of the tutorials to be tested
buildTutorialList(){
    nTutorial=0
    for d in $(find $COBRATutorialsPath -maxdepth 7 -type d)
    do
        if [[ "${d}" == *additionalTutorials* ]]; then
            continue  # if not a directory, skip
        fi

        # check for MLX files.
        for tutorial in ${d}/*.mlx
        do
            if ! [[ -f "$tutorial" ]]; then
                break
            fi
            let "nTutorial+=1"
            tutorials[$nTutorial]="$tutorial"
            echo " - ${tutorials[$nTutorial]}"
        done
    done
}

buildTutorialList

# preliminary tutorials list
declare -a tutorials=("tutorial_IO", "tutorial_COBRAconcepts", "tutorial_browseNetwork")

longest=0
for word in "${tutorials[@]}"
do
    len=${#word}
    if (( len > longest ))
    then
        longest=$len
    fi
done

header=`printf "%-${longest}s    %6s    %6s    %7s\n"  "Name" "passed" "failed" "time(s)"`
report="Tutorial report\n\n"
report+="$header\n"
report+=`printf '=%.0s' $(seq 1 ${#header});`"\n"
failure=0

# set time format to seconds
TIMEFORMAT=%R

nTutorial=0
nPassed=0

# loop through the tutorials
for tutorial in "${tutorials[@]}"
do
    tutorialDir=${tutorial%/*}
    tutorialName=${tutorial##*/}
    tutorialName="${tutorialName%.*}"

    msg="| Starting $tutorialName |"
    chrlen=${#msg}
    underline=`printf '=%.0s' $(seq 1 $chrlen);`
    echo "$underline"
    echo "$msg"
    echo "$underline"

    # time a process
    SECONDS=0;
    $ARTENOLIS_SOFT_PATH/MATLAB/$MATLAB_VER/bin/./matlab -nodesktop -nosplash -r "restoredefaultpath; addpath([pwd filesep 'tutorials' filesep '.artenolis']); runTutorial('$tutorialName'); delete(gcp);"
    CODE=$?
    procTime=$SECONDS

    msg="| Done executing $tutorialName! |"
    chrlen=${#msg}
    underline=`printf '=%.0s' $(seq 1 $chrlen);`
    echo "$underline"
    echo "$msg"
    echo "$underline"
    echo
    echo

    echo "Exit code: $CODE"

    if [ $CODE -ne 0 ]; then
        report+=`printf "%-${longest}s                x      %7.1f"  "$tutorial" "$procTime"`
    else
        report+=`printf "%-${longest}s     x                 %7.1f"  "$tutorial" "$procTime"`
        let "nPassed+=1"
    fi
    report+="\n"
    let "nTutorial+=1"
done

report+=`printf "\n  Passed:  %d/%d" "$nPassed" "$nTutorial"`
report+="\n\n"
printf "$report"

# exit
if [ $nPassed -ne $nTutorial ]; then
    exit 1
else
    exit $CODE
fi