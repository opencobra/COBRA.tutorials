## Continuous Integration of Tutorials:
### High level Overview
The whole idea of the continuous integration of the tutorials is that whenever a user contributes a tutorial in the format of a .mlx file on the [tutorials repo](https://github.com/opencobra/COBRA.tutorials) it should be converted to html and then rendered accordingly on the Cobratoolbox website in the tutorials section. This involves using MATLAB on a self-hosted server (King server) to generate the html file. This html is then pushed to the websites codebase repository which is the [./docs folder of the master branch](https://github.com/opencobra/cobratoolbox/tree/master/docs) of the main cobratoolbox repository.

GitHub actions is used to detect when a push (specifically .mlx push) is made to the tutorials repo. Then once the .mlx has been converted it is pushed to the gh-branch of the main repo. Again, GitHub actions can detect this push and configures the website to incorporate the extra tutorial. 

### Detailed Documentation
**Step 1: Pushing MLX files to the tutorials repository:**
To understand GitHub actions you need to look for the github workflow folder where you will find a .yml which contains all the details about the github action. The worflows can be found by navigating to ./.github/workflows/ . In the tutorials repo you will find a ‘[main.yml](https://github.com/opencobra/COBRA.tutorials/blob/master/.github/workflows/main.yml)’ file.

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
          repository: 'openCOBRA/COBRA.tutorials'
          token: ${{ secrets.GITHUB_TOKEN }}
          fetch-depth: 0
```


- Next, we have a series of ‘jobs’ to compute.
- The ‘runs-on’ parameter indicates where these jobs are computed. Here I specify it runs on ‘self-hosted’ because we need Matlab on King to run the .mlx to html. Generally, I would avoid using a self-hosted server but since Matlab is not an opensource programming language it needs to be ran a computer which has Matlab installed with a license.
- There are several steps to do in the jobs section. Here the first step is to checkout the source repo i.e. get all the details about the repo and the pushes made to the repo.

```
- name: Get All Changed Files
  id: files
  uses: jitterbit/get-changed-files@v1

- name: Find changed files
  id: getfile
  run: |
      files=$(echo "${{ steps.files.outputs.all }}" | tr ' ' '\n' | grep -E '\.mlx$')
      echo "::set-output name=file::$files"
```


- Here we have two more steps. The first step in this picture is used to find all the files that have been changed based on the most recent push.
- The next step is then used to find all the .mlx files that were pushed to the repository.

```
  - name: Give execute permission to the script
    run: chmod +x build.sh
    
  - name: Give execute permission to the other script  
    run: chmod +x setup.sh
```


The chmod command just makes the .sh files executable.

**Quick note on setup.sh and build.sh:**

• The setup.sh script automates the process of synchronizing .mlx files to the ghpages branch of the cobratoolbox GitHub repository. It requires three inputs: the repository identifier in owner/name format, a token for authentication, and the file path of the .mlx files to be synchronized. Upon execution, the script clones the cobratoolbox repository, configures git for automated operations, and targets aspecific directory within stable/tutorials/ to update. It clears this directory and copies the new .mlx files into it, ensuring that any changes are committed and pushed. This operation keeps the gh-pages branch of the cobratoolbox repository consistently updated with the latest .mlx files for documentation or tutorials.

• The build.sh script is designed for converting .mlx files to .html format and synchronizing them with the gh-pages branch of the cobratoolbox repository,. It takes three arguments: the repository identifier, a token for authentication, and the path of the .mlx file to be converted. Initially, the script converts the .mlx file to .html using MATLAB commands, assuming MATLAB is installed and accessible in the PATH. It then clones the target repository, sets up git with predefined user details, and switches to the gh-pages branch. The script creates a target directory within stable/tutorials/, copies the converted .html file into this directory, and finalizes by committing and pushing the changes.

• Both files can be found on the tutorial’s repository. Here are the links to [setup.sh](https://github.com/opencobra/COBRA.tutorials/blob/master/setup.sh) and [build.sh](https://github.com/opencobra/COBRA.tutorials/blob/master/build.sh)

```
  - name: Sync with Destination Repo
    run: |
      counter=0
      for file in ${{ steps.getfile.outputs.file }}; do
        if [ $counter -eq 0 ]
        then
          # This is the first iteration, handle the file differently
          ./setup.sh opencobra/cobratoolbox ${{ secrets.DEST_REPO_TOKEN }} $file 
          ./build.sh opencobra/cobratoolbox ${{ secrets.DEST_REPO_TOKEN }} $file 
        else
          ./build.sh openCOBRA/cobratoolbox ${{ secrets.DEST_REPO_TOKEN }} 
        fi
        counter=$((counter+1))
      done
    if: steps.getfile.outputs.file != ''
```

Here is the code to run the setup.sh and build.sh. We loop through all the .mlx files that were pushed. If it is the first file we are looking at we also run setup.sh to create the folder locations in the cobratoolbox – ghpages branch repository. Then afterwards build,sh is ran to convert the file to html and push to the created folder location

### Configuring the King Server

Go to this page of the repo to create a new self-hosted runner:

![image](https://github.com/opencobra/cobratoolbox/assets/68754265/05535af0-9ccf-4c38-9e79-512f738cc0f0)


By pressing the green new runner button, you are given easy instructions on how to set it up. You should have access to a terminal on King for this. To run the self-hosted runner nagivate to the folder you created it in and run ./run.sh to run the self-hosted runner.

You also need to make sure you have Matlab downloaded and working on the king server also. In the ‘[build.sh](https://github.com/opencobra/COBRA.tutorials/blob/master/build.sh)’ file the location of matlab is currently in my directory but you can add Matlab to another location and change the link to the location in the build.sh file.
