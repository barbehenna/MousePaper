#!/bin/bash

echo "This is only tested on a mac"
echo "This script will install things so that you can run the mouse catching app"

BREW=`which brew`
if [ -z "${BREW}" ]; then
	echo "homebrew not installed. \nInstalling homebrew...\n"
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
fi

R=`which R`
r=`which r`
if [ -z "${r}" -o -z "${R}" ]; then
	echo "R not installed. \nInstalling R...\n"
	brew tap homebrew/science
	brew install r
fi

echo "Installing R packages.\n"
R -e 'if(!require(shiny)){install.packages("shiny", repos="http://cran.us.r-project.org")}'
R -e 'if(!require(raster)){install.packages("raster", repos="http://cran.us.r-project.org")}'
R -e 'if(!require(Rcpp)){install.packages("Rcpp", repos="http://cran.us.r-project.org")}'

echo "To run the app, enter the command 'Rscript app.R' then paste the link shown into you browser."
