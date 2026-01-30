# Homework 01 - Data Input/Output with Python
Version: AA 2023-2024

## Introduction
This set of Python Notebooks and exercises are meant to be help in getting more familiar with the
Python functions used to read and manipulate datasets in two main formats used in multimessenger
science:

* CSV (Comma-separated values): ASCII-based text files
* FITS (Flexible Image Transport Sysyem), the standard in astrophysics

Later on, in future tutorials, we will see more formats, but for now let's focus on these two.

## What do you have to do?
You will find 2 sets of tutorial+exercises:
* **Set 1** - Reading and manipulating CSV files
* **Set 2** - Reading and manipulating FITS files

For each set, first look at the **tutorial**, check the code, run it and do modifications in order to be sure that you've understand how it works.

Then you can move to the **exercises**, where you will be asked to perform some tasks and solve some
real-life cases of data analysis.

## Structure of the package
The package contains the following directories

* *tutorials*: Contains Python Juyter notebooks for the tutorials
* *exercises*: Contains Python Juyter notebooks for the exercises
* *data*: Contains data required to run the tutorial and the exercises

## Getting the package
In order to get the package, you should run a git clone.
For instance 
```
git clone git@github.com:packageaddress
```
Where *git@github.com:packageaddress* has to be substituted with the proper link.

For instance, if your assigned package is called
*homework-01-data-input-output-with-python-myusername*, you should run:
```
git clone git@github.com:mmphyslab-pi/homework-01-data-input-output-with-python-myusername
```
(of course, you should substitute username with your Github username)
 
You can also look at the **Clone or Download** green button in the Github page of the package

**IMPORTANT** You might be asked for a password (if putting a HTTPS address), or you will be denied access (if using the ssh address). In order to create a SSH key and add it to Github, you can follow [these instructions](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account). Please follow **these steps** once for all so that you will not have to put the password every time you do a commit or a push.

## Working with the package
Of course, since this is a git repo, you can do what do you want (commit, push, create branches, etc).
However, if you want to modify and play with the tutorial code and not change the original, you can create a branch:
```
git branch dataio-mybranch
git checkout dataio-mybranch
git branch 
```
*Last command is to check if the branch is there...*


## Summarizing...
1) Set up the SSH passwordless access to Github;
2) Clone the repository;
3) Read the tutorials and play with them;
4) Do not forget to commit and push your changes often, in order to avoid losing your work in case of unexpected problems;
5) Work on the exercises;
6) When you have done the exercises, go to the final steps...


## At the end?
When you have completed your taks, please remember to commit and push adding a meaningful comment.
```
git commit -a -m "My name: Task Finished!" 
git push origin master dataio-mybranch
```


## Final thoughts
Some of you know Python well, some reasonably well, and some do not.
I hope these will be useful at different level for each of you to get more pratice on slightly advanced Python data analysis. 

These exercises will be useful for the various experiences of the course, so take advantage it it and enjoy!