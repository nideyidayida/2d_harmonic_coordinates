[![Build Status](https://travis-ci.org/NYUGP17/Assignment_5.svg?branch=master)](https://travis-ci.org/NYUGP17/Assignment_5)
# Assignment 5: Shape Deformation

This repository contains the viewer/painting code and data files you'll need for
assignment 5.

## Getting Started
To begin, clone the repository:
```
git clone https://github.com/NYUGP17/Assignment_5
```

Next, please refer to the [General Rules and Instructions](https://github.com/danielepanozzo/gp/raw/master/guidelines.pdf)
handout for instructions on installing LIBIGL and its dependencies.

## Building and Completing the Assignment
Once LIBIGL is set up (and pointed to by environment variable $LIBIGL_ROOT) you
should be able to build the viewer code:
```
mkdir build && cd build && cmake ..
make
```
Please report any problems you run into on this repository's Issues tab on
GitHub.

When the build completes successfully, begin implementing the missing blocks in
src/main.cpp as described by the assignment PDF.

## Submitting
When you finish the assignment, you will submit it by pushing it to a new
repository on our NYUGP17 organization.

1. Create a **private** repository in https://github.com/NYUGP17/ called
   **Assignment5_USER**, where USER is your github username that you entered in
   the survey.
2. Push your code to the repository:
```
git push https://github.com/NYUGP17/Assignment5_USER
```

## Travis-CI
Every submission must build on Linux before it can be graded/considered
complete. To check this, you will use Travis-CI, a tool for automatically
rebuilding your code each time you push it to GitHub.

We've already configured Travis-CI for this repository: Travis-CI works by
running the script in '.travis.yml' each time new commits are pushed. You will
need to follow the [getting started
instructions](https://travis-ci.com/getting_started) to sign into Travis-CI with
your GitHub account, grant it permission, and enable builds on your private
assignment repositories. Finally, you need to change the URL at the top of
this README.md file to point to your repository's status.
