# Experience - Imaging and Photometry in Optical with CCD

## Introduction
This experience is focused on the analysis of optical data taken with CCD and optical telescope, performing Image calibration and photometry.

As usual, you can find more information on the Official pages of these packages:

* [Official Page of the Numerical Python (numpy)](https://numpy.org/)
* [Official Page of the Astropy project](http://www.astropy.org)
* [Official Page of CCDProc](https://ccdproc.readthedocs.io/en/latest/)
* [Official Page of PhotUtils](https://photutils.readthedocs.io/en/stable/index.html)


As usual, here you can find a tutorial that is designed to prepare for the main analysis of the experience.

## What do you have to do?
First of all, have a look at the tutorial
* **Tutorial 01** - [Importing classes](tutorials/tutorial01-importing-classes.ipynb)

Once you have read these **tutorial**, check the code, run it and do modifications in order to be sure that you've understand how it works.

Then you can move to the description of the experience that you will find on Moodle, where you will be asked to perform some observations and data analysis on the Milky Way radio emission.

## Structure of the package
The package contains the following directories

* *tutorials*: Contains Python Juyter notebooks for the tutorials
* *data*: Contains data required to run the tutorials

## Getting the package
In order to get the package, you should run a git clone.
For instance 
```
git clone git@github.com:packageaddress
```
Where *git@github.com:packageaddress* has to be substituted with the proper link, as you have learned in previous exercises

You can also look at the **Clone or Download** green button in the Github page of the repository

**IMPORTANT** You might be asked for a password (if putting a HTTPS address), or you will be denied access (if using the ssh address). In order to create a SSH key and add it to Github, there are similar steps to those you have done already for gitlab. You can find more informtion at this [link on Github.com](https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh). Please follow **this step** once for all so that you will not have to put the password every time you do a commit or a push.

## Working with the package
Of course, since this is a git repo, you can do what do you want (commit, push, create branches, etc).
However, if you want to modify and play with the tutorial code and not change the original, you can create a branch:
```
git branch dataio-mybranch
git checkout dataio-mybranch
git branch 
```
Last command is to check if the branch is there...

**NOTE** Since this exercise is meant to be performed individually, you can avoid using branches and just push to the MASTER.


## Summarizing...
1) Set up the SSH passwordless access to Github;
2) Clone the repository;
3) Read the tutorials and play with them;
4) Do not forget to commit and push your changes often, in order to avoid losing your work in case of unexpected problems;
5) Work on the experience;
6) You should write a group report on it, you can save in the report directory
7) When you have done the exercise, go to the final steps...


## At the end?
When you have completed your taks, please remember to commit and push adding a meaningful comment.
```
git commit -a -m "My name: Task Finished!" 
git push origin master dataio-mybranch
```
