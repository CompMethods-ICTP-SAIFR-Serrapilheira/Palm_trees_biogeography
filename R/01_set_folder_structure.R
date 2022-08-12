####### Writing began in August 2022 by Thales
####### Written during the course Introduction to Scientific Computation
####### Written as part of the final evaluation of this same course
####### Serrapilheira ICTP/SAIRF QBIO program

####### This script sets the basic project folder structure

# Creating the data folder and the raw and processed data folders
# It is most important that the raw data remains untouched
# All the editions in data must be saved in the processed folder and be commented accordingly
if(!dir.exists("./data")){
  dir.create("./data")
}

if(!dir.exists("./data/raw")){
  dir.create("./data/raw")
}

if(!dir.exists("./data/processed")){
  dir.create("./data/processed")
}

# Creating the docs folder
# Here is where R Markdown (.Rmd) and .html knited files
if(!dir.exists("./docs")){
  dir.create("./docs")
}

# Creating the figs folder
# Here is where .png files will be saved
if(!dir.exists("./figs")){
  dir.create("./figs")
}

# Creating the R folder
# Here is where R scripts will be saved (but not those that are just defining functions)
if(!dir.exists("./R")){
  dir.create("./R")
}

# Creating the R folder
# Here is where R functions will be saved
if(!dir.exists("./fun")){
  dir.create("./fun")
}

# Creating the biblio folder
# Here is where .bib files with the references will be saved
if(!dir.exists("./biblio")){
  dir.create("./biblio")
}

# Creating the output folder
# Here is where outputs that are not figures or text documents will be saved
if(!dir.exists("./output")){
  dir.create("./output")
}

# Other folders maybe be created during the project, if this happen, just document it properly

# Creating a README file
usethis::use_readme_rmd()
