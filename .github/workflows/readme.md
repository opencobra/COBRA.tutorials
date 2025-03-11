## Continuous Integration of Tutorials:
### High level Overview
The whole idea of the continuous integration of the tutorials is that whenever a user contributes a tutorial in the format of a .mlx file on the [tutorials repo](https://github.com/opencobra/COBRA.tutorials) it should be converted to html and then rendered accordingly on the Cobratoolbox website in the tutorials section. This involves using MATLAB on a self-hosted server (King server) to generate the html file. This html is then pushed to the website's codebase repository which is the [./stable](https://github.com/opencobra/cobratoolbox/tree/gh-pages/stable) folder of the gh-pages branch of the main cobratoolbox repository.

GitHub actions is used to detect when a push (specifically .mlx push) is made to the tutorials repo. Then once the .mlx has been converted it is pushed to the stable folder of the gh-pages branch. Again, GitHub actions can detect this push and configures the website to incorporate the extra tutorial. 

### Detailed Documentation
**Step 1: Pushing MLX files to the tutorials repository:**
To understand GitHub actions you need to look for the github workflow folder where you will find a .yml which contains all the details about the github action. The worflows can be found by navigating to ./.github/workflows/. In the tutorials repo you will find a ‘[main.yml](https://github.com/opencobra/COBRA.tutorials/blob/master/.github/workflows/main.yml)’ file.

**What does main.yml do?**
Here is an explanation of each section of the .yml file. Pictures of the sections are added and an explanation is given beneath the picture.

```
on:
  push:
    branches: [ master ]
    paths:
    - '**.mlx'
```


This section of code basically means it will only run when a push is made to the master branch and one of the file types is a .mlx file. If not .mlx files are pushed, we don’t continue.

```
jobs:
  copy-changes:
    runs-on: self-hosted
    steps:
      - name: Checkout Source Repo
        uses: actions/checkout@v2
        with:
          repository: '${{ github.repository_owner }}/COBRA.tutorials'
          token: ${{ secrets.GITHUB_TOKEN }}
          fetch-depth: 0
```


- Next, we have a series of ‘jobs’ to compute.
- The ‘runs-on’ parameter indicates where these jobs are computed. Here, it runs on ‘self-hosted’ because we need Matlab on King to run the .mlx to html. Generally, we should avoid using a self-hosted server but since Matlab is not an opensource programming language it needs to be ran a computer which has Matlab installed with a license.
- There are several steps to do in the jobs section. Here the first step is to checkout the source repo i.e. get all the details about the repo and the pushes made to the repo.
  
```
- name: Get the repository's owner name
  run: |
    echo "REPO_OWNER=${{ github.repository_owner }}" >> $GITHUB_ENV
```
- Here we are getting the repository's owner name and storing it in the variable, REPO_OWNER. This variable will be used in further jobs.

```
- name: Clone the destination repository
  run: |
    rm -rf cobratoolbox
    echo "Cloning the destination repository: git@github.com:opencobra/cobratoolbox.git"
    git clone --depth 1 --branch gh-pages ssh://git@github.com/$REPO_OWNER/cobratoolbox.git
```
- Here cobratoolbox repo (only gh-pages branch) is cloned to push the generated .html pages in further steps 

```
- name: Get Changed mlx Files and sinc with destination repository
  id: getFile
  run: |
    changed_files=$(git diff --name-only HEAD~1 HEAD | grep '\.mlx' | tr '\n' ' ')
    # Setup virtual frame buffer
    for file in $changed_files; do
      if [[ $file != "" ]]; then
        echo "Processing: $file"
        ABSOLUTE_FILE_PATH=$(realpath "$file")
        HTML_FILE_PATH=$(echo "$ABSOLUTE_FILE_PATH" | sed 's/.mlx/.html/g')
        PDF_FILE_PATH=$(echo "$ABSOLUTE_FILE_PATH" | sed 's/.mlx/.pdf/g')
        M_FILE_PATH=$(echo "$ABSOLUTE_FILE_PATH" | sed 's/.mlx/.m/g')
        /usr/local/MATLAB/R2024a/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$ABSOLUTE_FILE_PATH', '$HTML_FILE_PATH')"
        /usr/local/MATLAB/R2024a/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$ABSOLUTE_FILE_PATH', '$PDF_FILE_PATH')"
        /usr/local/MATLAB/R2024a/bin/matlab -batch "matlab.internal.liveeditor.openAndConvert('$ABSOLUTE_FILE_PATH', '$M_FILE_PATH')"
        cd cobratoolbox
        TARGET_DIR="stable/tutorials/$(dirname "$file")"
        mkdir -p "$TARGET_DIR"
        echo "Copying the HTML, PDF, mlx and .m files to the target directory..."
        cp "$HTML_FILE_PATH" "$TARGET_DIR/"
        cd ../
      fi
    done
```
- Here we actually do the conversions and copy the .HTML file to cobratoolbox repository. The fourth line in this picture (changed_files=$(git diff --name-only HEAD~1 HEAD | grep '\.mlx' | tr '\n' ' ')) is used to find all the .mlx files that have been changed based on the most recent push.
- Then for each modified .mlx file, an HTML file, a pdf file, and a Matlab code file are created and stored in the same directory where the .mlx file is located. We are using MATLAB R2024a here, if the king server has some other version, then this step has to be modified accordingly. Further, HTML file alone is copied to the cobratoolbox repo.


```
- name: Pushing the changes to both COBRA.tutorials and cobratoolbox repos
  run: |
    # Set up git config
    echo "Setting up git config..."
    git config --global user.name "github-actions[bot]"
    git config --global user.email "github-actions[bot]@users.noreply.github.com"
    # Add, commit, and push the changes
    cd cobratoolbox
    git pull
    git add .
    git commit -m "Sync files from source repo" || echo "No changes to commit"
    git push
    cd ..
    rm -rf cobratoolbox
    git add .
    git commit -m "created .HTML, .pdf and .m files"
    git push origin master
    echo "Script execution completed."
```

The converted files are further pushed to the cobratoolbox and COBRA.tutorial repo. Note that only .html pages are pushed to cobratoolbox repository and all the other formats are stored in COBRA.tutorial repo

### Configuring the King Server

Go to this page of the repo to create a new self-hosted runner:

![image](https://github.com/opencobra/cobratoolbox/assets/68754265/05535af0-9ccf-4c38-9e79-512f738cc0f0)


By pressing the green new runner button, we are given easy instructions on how to set it up. We should have access to a terminal on King for this. To run the self-hosted runner nagivate to the folder we created it in and run ./run.sh to run the self-hosted runner.

We also need to make sure you have Matlab downloaded and working on the king server also.
