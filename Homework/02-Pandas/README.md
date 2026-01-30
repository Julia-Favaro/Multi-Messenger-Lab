# Homework 02 - Manipulation of Astrophysical data with Astropy and Pandas

Version: AA 2023-2024

## Introduction
This set of Python Notebooks and exercises are meant to be help in getting more familiar with the
Python libraries Numpy, Astropy and Pandas, that are used to manipulate tables and data in the astrophysical contect.
In particular you can find more information on the Official pages of these packages:

* [Official Page of the Numerical Python (numpy)](https://numpy.org/)
* [Official Page of the Astropy project](http://www.astropy.org)
* [Official Page of Pandas](https://pandas.pydata.org)

After having looked at a set of tutorials you will develop a small project to look for spatial correlations among source catalogs.

## What do you have to do?
You will find one sets of tutorials and a final exercise:
* **Tutorial 01** - [Basics of Numerical Python](tutorials/tutorial01-numpy-basics.ipynb)
* **Tutorial 02** - [Managing times with Astropy](tutorials/tutorial02-astropy-times.ipynb)
* **Tutorial 03** - [Managing coordinates with Astropy](tutorials/tutorial03-astropy-coordinates.ipynb)
* **Tutorial 04** - [Basics of Pandas](tutorials/tutorial04-pandas.ipynb)
* **Tutorial 05** - [Putting all together: Pandas and FITS files](tutorials/tutorial05-pandas-fits.ipynb)

Once you have read these **tutorial**, check the code, run it and do modifications in order to be sure that you've understand how it works.

Then you can move to the [**exercise**](exercises/exercise01-x-gamma-catalogs.ipynb), where you will be asked to perform some tasks and solve a real-life cases of data analysis, in this case looking for spatial correlations among X-ray and gamma-ray sources.

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
git clone git@github.com:mmphyslab-pi/homework-02-python-astrodata-template-myusername
```
(of course, you should substitute username with your Github username)
 
You can also look at the **Clone or Download** green button in the Github page of the package

**IMPORTANT** You might be asked for a password (if putting a HTTPS address), or you will be denied access (if using the ssh address).  In order to create a SSH key and add it to Github, you can follow [these instructions](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
Please follow **these steps** once for all so that you will not have to put the password every time you do a commit or a push.

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
5) Work on the exercise;
6) When you have done the exercise, go to the final steps...


## At the end?
When you have completed your taks, please remember to commit and push adding a meaningful comment.
```
git commit -a -m "My name: Task Finished!" 
git push origin master dataio-mybranch
```


## Final thoughts
Some of you know Python well, some reasonably well, and some do not.
As for the previous homework, I hope this will be useful at different level for each of you to get more pratice on slightly advanced Python data analysis. 

These exercises will be useful for the various experiences of the course, so take advantage it it and enjoy!