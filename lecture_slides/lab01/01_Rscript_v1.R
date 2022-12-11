


# This is a comment
# We will paste in the example R script from:
# http://phylo.wikidot.com/introduction-to-r-pcms

# Instructions:
# 1. Open a *plain-text* editor (Mac: TextWrangler, BBedit, 
# R.app, RStudio) (Windows: Notetab++, RStudio)
# 
# 2. Copy/paste in the example script
# 
# 3. Save it as a .R file, in a 
#    directory like REU_example
# 
# 4. Have R open on the left, and the 
# text file open on the right

##############################################################
# =============================================
# Introduction to R and Phylogenetic Comparative Methods
# =============================================
# by Nick Matzke (and whoever else adds to this PhyloWiki page)
# Copyright 2014-infinity
# matzkeATberkeley.edu
# Last update: May 2019
#
# Please link/cite if you use this, email me if you have 
#   thoughts/improvements/corrections.
#
##############################################################
#
# Reference: Matzke, Nicholas J. (2020). "Introduction to R and Phylogenetic Comparative Methods." 
# Lab exercise for Transmitting Science, Barcelona, Spain.
# December 7, 2020
# 
#
# Previous versions:
# 
# Reference: Matzke, Nicholas J. (2019). "Introduction to R and Phylogenetic Comparative Methods." 
# Lab exercise for Bioinformatics 702 (BIOINF702) at the University of Auckland.
# May 7, 2019.
# 
# Free to use/redistribute under:
# Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0) 
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the above license, linked here:
# 
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# 
# Summary:
#
# You are free:
#
#   * to Share — to copy, distribute and transmit the work
#   * to Remix — to adapt the work
#
# Under the following conditions:
#
#   * Attribution — You must attribute the work in the manner 
#     specified by the author or licensor (but not in any way that 
#     suggests that they endorse you or your use of the work).
#   * Noncommercial — You may not use this work for commercial purposes. 
#
#   * Share Alike — If you alter, transform, or build upon this work,
#     you may distribute the resulting work only under the same or
#     similar license to this one. 
#
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# 
###################################################################

#######################################################
# CHAPTER 1: PROPAGANDA FOR R
# 
# R is a programming language designed primarily for 
# data analysis and statistics.
#
# The big advantages of R are:
#
# 1. It is free.
# 2. It is easy.
#
# Point #2 sometimes takes some convincing, especially
# if you haven't programmed before. But, trust me, R
# is WAY easier than ANY other programming language
# I have ever tried, which you could also do serious
# science with.
#
# MATLAB is probably the only other competitor for ease 
# of use and scientific ability, but Matlab costs 
# hundreds of dollars, and hundreds of dollars more for
# the various extensions (for e.g. statistics, image 
# analysis, etc.).  This works great when your institution
# has a site license for Matlab, but it suck when you 
# move to a new school/job.
#
# R is easy because most of the "computer science 
# details" -- how to represent numbers and other 
# objects in the computer as binary bits/bytes,
# how to manage memory, how to cast and type variables,
# blah blah blah, are done automatically behind the 
# scenes. 
# 
# This means almost anyone can get going with R in 
# minutes, by just typing in commands and not having
# to spend days learning the difference between a 
# short and long integer, blah blah blah.
#
# That said, the cost of this automation is that R
# is slower than other programming languages.  However, 
# this doesn't matter for common, basic sorts of 
# statistical analyses -- say, linear regression with 
# 1,000 data observations.  It DOES matter if you are 
# dealing with huge datasets -- say, large satellite
# images, or whole genomes.
#
# In these situations, you should use specialist 
# software, which is typically written in Python
# (for manipulating textual data, e.g. genome files)
# or Java, C, or C++ (for high-powered computing).
#
# (Although, in many situations, the slow parts of 
#  R can be re-programmed in C++, and accessed from
#  R.)
#
# R is also pretty bad for large, complex programming
# projects.  Python and C++ are "object-oriented."
# In computer-programming, "objects" help organize
# your data and tasks.  For example, if you are 
# writing a video game, you might want to program 
# many different monsters. However, you don't want to
# re-program the behavior of each monster from scratch.
# Instead, you create a general object, "monster", and 
# give it attributes (speed, armor, etc.).  The "monster"
# object takes inputs (like what enemies are close to 
# it) and produces outputs (motion or attacks in a 
# certain direction).
#
# Each specific type of monster would be an instance 
# of the monster class of objects.  Each individual
# monster of a specific type would be its own object,
# keeping track of hit points, etc.
#
# You can see that, for serious programming, this 
# object-oriented style would be the way to go. Therefore,
# "real" computer-science classes teach you this way
# of programming. This is great if you want to 
# go work in the video game industry and devote your
# life to coding.
#
# However, if you just want to plot some data and 
# run some statistical tests and do some science, 
# you don't want to have to go through a bunch of 
# rigamarole first. You just want to load the data
# and plot it and be done.  This is what R is for.
#
#######################################################

#######################################################
# CHAPTER 2: GETTING R
#######################################################
# 
# R is free and available for all platforms. You can 
# download it here.:
# 
# http://www.r-project.org/
# 
# Tip for free, scientific software: 
# 
# Unless you are doing something expert, you will want 
# the "binary" file rather than the source code.
#
# Programmers write source code in text files.
# 
# A compiler program turns this into a "binary" which 
# actually executes (runs) on a computer.
# 
# Compiling from source code can take minutes or hours, 
# and sometimes will crash if your computer & compiler
# are not set up right.
# 
# A binary should just work, once you have installed it,
# assuming you've got the binary for your machine.
#
# ASSIGNMENT: Once you have R installed (it appear in 
# "Applications" on a Mac, or "Program Files" on a 
# Windows machine), open it to make sure it works. 
# Then, return to this tutorial.
# 
########################################################

#######################################################
# CHAPTER 3: GET A *PLAIN*-TEXT EDITOR
#######################################################
#
# Many people make the mistake of typing commands 
# into R, but not saving those commands.
# 
# *ALWAYS* SAVE YOUR COMMANDS IN A TEXT FILE!!
# *ALWAYS* SAVE YOUR COMMANDS IN A TEXT FILE!!
# *ALWAYS* SAVE YOUR COMMANDS IN A TEXT FILE!!
#
# Got it?  Good.
#
# The next mistake people make is to use Word or 
# some other monstrosity to save their commands. 
# You can do this if you want, but the formatting 
# etc. just gets in the way.  
#
# Find or download a PLAIN-TEXT editor (aka ASCII
# text editor). Common examples:
#
# Mac: TextWrangler (free) or BBedit 
#
# Windows: Notepad (free, search Programs) or Notetab
# 
# Or: versions of R that have a GUI (GUI=Graphical User
# Interface) also have a built-in editor.
# 
#
# WHY SAVE YOUR COMMANDS?
# 
# Because you can come back in 6 months and run the 
# same analysis again, just by pasting the commands
# back into R.
#
# Trust me, this is MUCH better than trying to remember
# what buttons to click in some software. 
# 
# And, anytime
# you need to do something more than a few times, 
# it gets super-annoying to click all of the buttons
# again and again. 
# 
# This is why most serious scientific software is 
# command-line, rather than menu-driven.
# 
#
# HOW TO TAKE NOTES IN R SCRIPTS
# 
# Put a "#" symbol in front of your comments. Like I 
# did here. COMMENTS ARE GOOD! COMMENT EVERYTHING!
# 
#
# ASSIGNMENT: Once you've found a plain-text editor, 
# return to this tutorial.
#######################################################

#######################################################
# CHAPTER 4: R BASICS
#######################################################
# 
# There are two major hurdles in learning R:
# 
# 1. Getting/setting your working directory.
#
# 2. Loading your data
#
# 3. Learning the commands to do what you want.
#
# Points #1 and #2 are easy to learn -- just don't
# forget!  You can never get anything significant
# done in R if you can't get your data loaded.
# 
# Point #3 -- No one knows "all" of R's commands. As
# we see, every package and function creates 
# additional commands.
# 
# Your goal is just to learn the basics, and then learn
# how to find the commands you need.
# 
# ASSIGNMENT: Type/paste in each of the commands below
# into your text file, then into R. Take notes as 
# you go.
#######################################################

#######################################################
# Working directories:
#
# One of the first things you want to do, usually, is 
# decide on your working directory.
# 
# You should create a new directory using:
#
# Mac: Finder
# Windows: Windows Explorer (or File Manager or
# whatever it's called these days)
#
# ROOLZ FOR FILES AND DIRECTORIES IN R
# 
# 1. Put the directory somewhere you will find it 
#    again.
#
# 2. Never use spaces in filenames.
# 
# 3. Never use spaces in directory names.
# 
# 4. Never use spaces in anything involving 
#    files/directories.
#
# 5. Never!  It just causes problems later. The 
#    problems are fixable, but it's easier to 
#    just never use spaces.  
# 
# 6. Use underscore ("_") instead of spaces.
# 
#
# FINDING MAC/WINDOWS DIRECTORIES IN R
#
# Usually, you can drag/drop the file or directory
# into R to see the full path to the file.
# Copy this into the 'here' in wd="here", below.
#
# 
# CHANGE FILE SETTINGS IN MACS/WINDOWS
# 
# Modern Macs/Windows hide a lot of information 
# from you.  This makes life easier for John Q. Public,
# but makes it harder for scientists.  
# 
# Good preferences for your file viewer:
#
# * Turn ON viewing file extensions (.txt, .docx, etc.)
# * Turn ON viewing of hidden files 
# * Change file viewing to "list" format
#
# See Preferences in Mac Finder or Windows Explorer.
# 
########################################################

# On my Mac, this is a working directory I have chosen
# (change it to be yours)
wd = "/Users/nickm/Desktop/Rintro/"
wd = "~/Desktop/Rintro/"
wd="/drives/GDrive/REU_example"

# On an older PC, you might have to specify paths like this:
#wd = "c:\\Users\\nick\\Desktop\\_ib200a\\ib200b_sp2011\\lab03"

# setwd: set working directory
setwd(wd)

# getwd: get working directory
getwd()

# list.file: list the files in the directory
list.files()

#######################################################
# PLAYING WITH R
#######################################################

# (Preliminary: this might be useful; uncomment if so)
# options(stringsAsFactors = FALSE)

# concatentate to a list with c()
student_names = c("Nick", "Hillary", "Sonal")

# describe what the function "c" does:
#
grade1 = c(37, 100, 60)
grade2 = c(43, 80, 70)
grade3 = c(100, 90, 100)
grade1
grade2
grade3
print(grade3)

# column bind (cbind)
temp_table = cbind(student_names, grade1, grade2, grade3)
class(temp_table)

# convert to data frame
grade_data = as.data.frame(temp_table)
class(grade_data)

# Don't convert to factors
grade_data = as.data.frame(temp_table, stringsAsFactors=FALSE)

# add column headings
col_headers = c("names", "test1", "test2", "test3")
names(grade_data) = col_headers
print(grade_data)

# change the column names back
old_names = c("student_names", "grade1", "grade2", "grade3")
names(grade_data) = old_names

grade_data$grade1

# Let's calculate some means
# mean of one column
mean(grade_data$grade1)

# R can be very annoying in certain situations, e.g. treating numbers as character data
# What does as.numeric do?
#
as.numeric(grade_data$grade1)
grade_data$grade1 = as.numeric(as.character(grade_data$grade1))
grade_data$grade2 = as.numeric(as.character(grade_data$grade2))
grade_data$grade3 = as.numeric(as.character(grade_data$grade3))
print(grade_data)

# mean of one column
mean(grade_data$grade1)

# apply the mean function over the rows, for just the numbers columns (2, 3, and 4)
apply(X=grade_data[ , 2:4], MARGIN=2, FUN=mean)

# Why doesn't this work?
mean(grade_data)
# What caused the warning message in mean(grade_data)?

# How about this?
colMeans(grade_data[,2:4])

# How about this?
colMeans(grade_data[,2:4])

# More functions
sum(grade_data$grade1)
median(grade_data$grade1)

# standard deviation
apply(X=grade_data[ , 2:4], MARGIN=1, FUN=sd)

# store st. dev and multiply by 2
mean_values = apply(grade_data[ , 2:4], 1, mean)
sd_values = apply(grade_data[ , 2:4], 1, sd)
2 * sd_values

# print to screen even within a function:
print(sd_values)

# row bind (rbind)
grade_data2 = rbind(grade_data, c("means", mean_values), c("stds", sd_values))

#######################################################
# GETTING DATA
#######################################################
#
# Let's download some data.  Francis Galton was one 
# of the founders of statistics.  He was also 
# the cousin of Charles Darwin. Galton invented the 
# term "regression".  These days, "regression" means
# fitting the best-fit line to a series of x and y
# data points.
# 
# But, why is the weird term "regression" used for this?
# What is regressing?
#
# Let's look at Galton's original dataset: the heights
# of parents and children.
#
# Use your web browser to navigate here:
#
# http://www.randomservices.org/random/data/Galton.html
# 
# ...and save "Galton's height data" as Galton.txt 
# (right-click, save) to your
# working directory.
#
# After doing this, double-click on Galton.txt and 
# view the file, just to see what's in there.
#
#######################################################

# Before proceeding, double-check that your data file
# is in the working directory:
getwd()
list.files()

# Let's store the filename in a variable
# 
# Note: In Nick's head:
# 
# "wd" means "working directory"
# "fn" means "filename"
# 
#wd = "/drives/Dropbox/_njm/__packages/Rintro/"
#setwd(wd)
fn = "Galton2.txt"

# Now, read the file into a data.frame
heights = read.table(file=fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Now, look at "heights"
heights

# Whoops, that went by fast!  Let's just look at the 
# top of the data table
head(heights)

# Let's get other information on the data.table

# Column names
names(heights)

# Dimensions (rows, columns)
dim(heights)

# Class (data.frame, matrix, character, numeric, list, etc.)
class(heights)

# The heights data is the adult height of a child (in inches),
# and the "midparent" height -- the mean of the two parents.

# QUESTION: Do the means of parent and child height differ?

# Means
colMeans(heights)

colMeans(heights[,c(-1,-4)])

# Standard deviations
apply(X=heights[,c(-1,-4)], MARGIN=2, FUN=sd)

# Min & Max
apply(X=heights[,c(-1,-4)], MARGIN=2, FUN=min)
apply(X=heights[,c(-1,-4)], MARGIN=2, FUN=max)

# They seem pretty close, but let's do a test

# Make sure numbers columns are numeric
heights$Family = as.numeric(heights$Family)
heights$Father = as.numeric(heights$Father)
heights$Height = as.numeric(heights$Height)
heights$Kids = as.numeric(heights$Kids)

# Let's add the Midparent column
heights[,c("Father","Mother")]

# Take the mean of Father and Mother columns, store in column "Midparent"
heights$Midparent = apply(X=heights[,c("Father","Mother")], MARGIN=1, FUN=mean)

# View the new column
head(heights)

# Population Mean Between Two Independent Samples
# http://www.r-tutor.com/elementary-statistics/inference-about-two-populations/population-mean-between-two-independent-samples

# (change "Child" to "Height")
ttest_result1 = t.test(x=heights$Midparent, y=heights$Height, paired=FALSE, alternative="two.sided")
ttest_result1

# But wait, this test assumes that the samples from each population 
# are independent. Do you think parent heights and child heights are 
# independent?

# Probably not.  Actually, these samples are paired, so let's
# check that:

# Population Mean Between Two Matched Samples
# http://www.r-tutor.com/elementary-statistics/inference-about-two-populations/population-mean-between-two-matched-samples

ttest_result2 = t.test(x=heights$Midparent, y=heights$Height, paired=TRUE, alternative="two.sided")
ttest_result2

# Compare the two:
ttest_result1
ttest_result2

# Interestingly, it looks like parents are slightly taller than the children!
# 
# Is this statistically significant?
#
# But is it a large effect?  Is it *practically* significant?
#

# Let's plot the histograms
hist(heights$Midparent)
hist(heights$Height)

# That's a little hard to compare, due to the different 
# automated scaling of the x-axis.

# Let's fix the x-axis to be (5 feet, 7 feet)
xlims = c(5*12, 7*12)
hist(heights$Midparent, xlim=xlims)
hist(heights$Height, xlim=xlims)

# And fix the y-axis
# Let's fix the y-axis to be (0, 220)
ylims = c(0, 220)
hist(heights$Midparent, xlim=xlims, ylim=ylims)
hist(heights$Height, xlim=xlims, ylim=ylims)

# Let's plot the means and 95% confidence intervals on top

# Midparent values
hist(heights$Midparent, xlim=xlims, ylim=ylims)

# Plot the mean
abline(v=mean(heights$Midparent), lty="dashed", lwd=2, col="blue")

# Plot the 95% confidence interval (2.5% - 97.5%)
CI_025 = mean(heights$Midparent) - 1.96*sd(heights$Midparent)
CI_975 = mean(heights$Midparent) + 1.96*sd(heights$Midparent)
abline(v=CI_025, lty="dotted", lwd=2, col="blue")
abline(v=CI_975, lty="dotted", lwd=2, col="blue")

# Child values
hist(heights$Height, xlim=xlims, ylim=ylims)

# Plot the mean
abline(v=mean(heights$Height), lty="dashed", lwd=2, col="blue")

# Plot the 95% confidence interval (2.5% - 97.5%)
CI_025 = mean(heights$Height) - 1.96*sd(heights$Height)
CI_975 = mean(heights$Height) + 1.96*sd(heights$Height)
abline(v=CI_025, lty="dotted", lwd=2, col="blue")
abline(v=CI_975, lty="dotted", lwd=2, col="blue")

# Let's put it all in a nice PDF format to save it

# Open a PDF for writing
pdffn = "Galton_height_histograms_v1.pdf"
pdf(file=pdffn, width=8, height=10)

# Do 2 subplots
par(mfrow=c(2,1))

# Midparent values
hist(heights$Midparent, xlim=xlims, ylim=ylims, xlab="height (inches)", ylab="Count", main="Midparent heights")

# Plot the mean
abline(v=mean(heights$Midparent), lty="dashed", lwd=2, col="blue")

# Plot the 95% confidence interval (2.5% - 97.5%)
CI_025 = mean(heights$Midparent) - 1.96*sd(heights$Midparent)
CI_975 = mean(heights$Midparent) + 1.96*sd(heights$Midparent)
abline(v=CI_025, lty="dotted", lwd=2, col="blue")
abline(v=CI_975, lty="dotted", lwd=2, col="blue")

# Child values
hist(heights$Height, xlim=xlims, ylim=ylims, xlab="height (inches)", ylab="Count", main="Child heights")

# Plot the mean
abline(v=mean(heights$Height), lty="dashed", lwd=2, col="blue")

# Plot the 95% confidence interval (2.5% - 97.5%)
CI_025 = mean(heights$Height) - 1.96*sd(heights$Height)
CI_975 = mean(heights$Height) + 1.96*sd(heights$Height)
abline(v=CI_025, lty="dotted", lwd=2, col="blue")
abline(v=CI_975, lty="dotted", lwd=2, col="blue")

# Close the PDF writing
dev.off()

# Write a system command as a text string
cmdstr = paste("open ", pdffn, sep="")
cmdstr

# Send the command to the computer system's Terminal/Command Line
system(cmdstr)
# The PDF should hopefully pop up, e.g. if you have the free Adobe Reader

# The difference in means is very small, even though it appears to be 
# statistically significant. 
#
# This is a VERY IMPORTANT lesson: 
#
# "statistically significant" DOES NOT ALWAYS MEAN "practically "significant",
# "interesting", "scientifically relevant", etc.
# 
# 
# The difference may have to do with:
# 
# * Galton's 'method' of dealing with the fact that 
# male and female children have different average heights -- 
# he multiplied the female heights by 1.08!
#
# * Different nutrition between the generations
#
# * Maybe the adult children weren't quite all fully grown
#
# * Chance rejection of the null
#
# Who knows?

# You may have noticed that the standard deviations look to be 
# a lot different.  Can we test for this?

# Yes! The null hypothesis is that the ratio of the 
# variances is 1:
Ftest_result = var.test(x=heights$Midparent, y=heights$Height, ratio=1, alternative="two.sided")
Ftest_result

# We get extremely significant rejection of the null.  What is 
# the likely cause of the lower variance in the midparent data?

#
# For the complex story of Galton's original data, see:
# 
# http://www.medicine.mcgill.ca/epidemiology/hanley/galton/
#
# James A. Hanley (2004). 'Transmuting' women into men: 
# Galton's family data on human stature. The American Statistician, 58(3) 237-243. 
# http://www.medicine.mcgill.ca/epidemiology/hanley/reprints/hanley_article_galton_data.pdf
# 
# BTW, Galton was both a genius, and promoted some deeply flawed ideas
# like eugenics:
# http://isteve.blogspot.com/2013/01/regression-toward-mean-and-francis.html
# 

# We noted before that child and parent heights might not be 
# independent.  Let's test this!

# QUESTION: is there a relationship?

# Start by plotting the data:
plot(x=heights$Midparent, y=heights$Height)

# It looks like there is a positive relationship: 
# taller parents have taller children.

# However, it's a little bit hard to tell for 
# sure, because Galton's data is only measured
# to the half-inch, so many dots are plotting
# on top of each other.  We can fix this by
# "jittering" the data:

# Plot the data, with a little jitter
plot(x=jitter(heights$Midparent), y=jitter(heights$Height))

# It looks like there's a positive relationship, which makes
# sense.  Can we confirm this with a statistical test?

# Let's build a linear model (lm)
lm_result = lm(formula=Height~Midparent, data=heights)
lm_result

# This just has the coefficients, this doesn't tell us much
# What's in the linear model? A list of items:
names(lm_result)

# See the statistical results
summary(lm_result)

# Analysis of variance (ANOVA)
anova(lm_result)

# You can get some standard diagnostic regression plots with:
plot(lm_result)

# Let's plot the regression line on top of the points
intercept_value = lm_result$coefficients["(Intercept)"]
slope_value = lm_result$coefficients["Midparent"]

# Plot the points
plot(x=jitter(heights$Midparent), y=jitter(heights$Height))

# Add the line
abline(a=intercept_value, b=slope_value, col="blue", lwd=2, lty="dashed")

# It's a little hard to tell if the slope is 1:1 or not,
# Because the x-axis and y-axis aren't the same
# Let's fix this

# Plot the points
xlims = c(5*12, 6.5*12)
ylims = c(5*12, 6.5*12)
plot(x=jitter(heights$Midparent, factor=3), y=jitter(heights$Height, factor=3), xlab="Midparent height", ylab="Child height", xlim=xlims, ylim=ylims)
title("Galton's height data")

# Add the regression line
abline(a=intercept_value, b=slope_value, col="blue", lwd=2, lty="dashed")

# Add the 1:1 line
abline(a=0, b=1, col="darkgreen", lwd=2, lty="dashed")

# Is the slope statistically different from 1:1?

# We can test this by subtracting a 1:1 relationship from the data, and seeing if 
# the result has a slope different from 0
child_minus_1to1 = heights$Height - (1/1*heights$Midparent)
heights2 = heights
heights2 = cbind(heights2, child_minus_1to1)

# Let's build a linear model (lm)
lm_result2 = lm(formula=child_minus_1to1~Midparent, data=heights2)
lm_result2

# This just has the coefficients, this doesn't tell us much
# What's in the linear model? A list of items:
names(lm_result2)

# See the statistical results
summary(lm_result2)

# Analysis of variance (ANOVA)
anova(lm_result2)

# You can get some standard diagnostic regression plots with:
plot(lm_result2)

# Let's plot the regression line on top of the points
intercept_value = lm_result2$coefficients["(Intercept)"]
slope_value = lm_result2$coefficients["Midparent"]

# Plot the points
plot(x=jitter(heights2$Midparent), y=jitter(heights2$child_minus_1to1), xlim=xlims, xlab="Midparent heights", ylab="Child heights minus 1:1 line", main="Relationship after subtracting 1:1 line")

# Add the regression line
abline(a=intercept_value, b=slope_value, col="blue", lwd=2, lty="dashed")

# Add the expected line if the relationship was 1:1
abline(a=0, b=0, col="darkgreen", lwd=2, lty="dashed")

# Yep, the relationship is definitely different than 1:1

# Why is the relationship between parent height and offspring
# height LESS THAN 1:1???
#
# Why do tall parents tend to produce offspring shorter
# than themselves?   Why does height seem to "regress"?
# What about the children of short parents?  Do they 
# 'regress'?
# 
# What are possible statistical consequences/hazards of this?
#
# Why is all of this rarely explained when regression 
# is taught?
# 

#######################################################
# CHAPTER 5: MAKE YOUR OWN FUNCTIONS, AND DO MAXIMUM LIKELIHOOD
#
# R has many good functions, but it is easy to make your 
# own!  In fact, this is necessary for some applications.
#
#######################################################

# Let's consider some coin-flip data.  
#
# Here are 100 coin flips:

coin_flips = c('H','T','H','T','H','H','T','H','H','H','T','H','H','T','T','T','T','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','T','T','T','H','T','T','T','H','T','T','T','H','H','H','T','T','H','H','H','T','H','H','H','T','T','H','H','H','H','H','H','H','T','T','H','H','H','H','T','T','H','H','H','T','T','H','H','H','H','H','H','T','T','T','H','H','H','H','H','H','T','H','T','H','H','T','T')

coin_flips

# What is your guess at "P_heads", the probability of heads?
#
# What do you think the Maximum Likelihood (ML) estimate would be?
#

# In the case of binomial data, we actually have a formula to calculate
# the ML estimate:

# Find the heads
heads_TF = (coin_flips == "H")
heads_TF

# Find the tails
tails_TF = (coin_flips == "T")
tails_TF

numHeads = sum(heads_TF)
numHeads

numTails = sum(tails_TF)
numTails

numTotal = length(coin_flips)
numTotal

# Here's the formula:
P_heads_ML_estimate = numHeads / numTotal
P_heads_ML_estimate

# Well, duh, that seems pretty obvious.  At least it would have been, if we 
# weren't thinking of coins, where we have a strong prior belief that the
# coin is probably fair.

# What does it mean to say that this is "maximum likelihood" estimate of P_heads?
# 
# "Likelihood", in statistics, means "the probability of the data under the model"
#
# This technical definition has some interesting consequences:
# 
# * Data has likelihood, models do not.
# * The likelihood of noise in my attic, under the model that grelims
#   are having a party up there, is 1.

# Let's calculate the probability of the coin flip data under the 
# hypothesis/model that P_heads is 0.5

# We'll be very inefficient, and use a for-loop, and
# if/else statements

# Loop through all 100 flips
# Make a list of the probability of 
# each datum
P_heads_guess = 0.5

# Empty list of probabilities
probs_list = rep(NA, times=length(coin_flips))
probs_list

for (i in 1:length(coin_flips))
    {
    # Print an update
    cat("\nAnalysing coin flip #", i, "/", length(coin_flips), sep="")

    # Get the current coin flip
    coin_flip = coin_flips[i]

    # If the coin flip is heads, give that datum
    # probability P_heads_guess.
    # If tails, give it (1-P_heads_guess)

    if (coin_flip == "H")
        {
        probs_list[i] = P_heads_guess
        } # End if heads

    if (coin_flip == "T")        
        {
        probs_list[i] = (1-P_heads_guess)
        } # End if tails
    } # End for-loop

# Look at the resulting probabilities
probs_list

# We get the probability of all the data by multiplying
# all the probabilities
likelihood_of_data_given_P_heads_guess1 = prod(probs_list)
likelihood_of_data_given_P_heads_guess1

# That's a pretty small number!  You'll see that it's 
# just 0.5^100:
0.5^100

# A probability of 0.5 is not small, but multiply it 
# 100 values of 0.5 together, and you get a small value.
# That's the probability of that specific sequence of 
# heads/tails, given the hypothesis that the true
# probability is P_heads_guess.

# Let's try another probability:

# Loop through all 100 flips
# Make a list of the probability of 
# each datum
P_heads_guess = 0.7

# Empty list of probabilities
probs_list = rep(NA, times=length(coin_flips))
probs_list

for (i in 1:length(coin_flips))
    {
    # Print an update
    cat("\nAnalysing coin flip #", i, "/", length(coin_flips), sep="")

    # Get the current coin flip
    coin_flip = coin_flips[i]

    # If the coin flip is heads, give that datum
    # probability P_heads_guess.
    # If tails, give it (1-P_heads_guess)

    if (coin_flip == "H")
        {
        probs_list[i] = P_heads_guess
        } # End if heads

    if (coin_flip == "T")        
        {
        probs_list[i] = (1-P_heads_guess)
        } # End if tails
    } # End for-loop

# Look at the resulting probabilities
probs_list

# We get the probability of all the data by multiplying
# all the probabilities
likelihood_of_data_given_P_heads_guess2 = prod(probs_list)
likelihood_of_data_given_P_heads_guess2

# We got a different likelihood. It's also very small.
# But that's not important. What's important is, 
# how many times higher is it?

likelihood_of_data_given_P_heads_guess2 / likelihood_of_data_given_P_heads_guess1

# Whoa!  That's a lot higher!  This means the coin flip data is 54 times more
# probable under the hypothesis that P_heads=0.7 than under the 
# hypothesis that P_heads=0.5.

# Maximum likelihood: You can see that the BEST explanation of the data 
# would be the one with the value of P_heads that maximized the probability 
# of the data.  This would be the Maximum Likelihood solution.

# We could keep copying and pasting code, but that seems annoying.  Let's make a function 
# instead:

# Function that calculates the probability of coin flip data
# given a value of P_heads_guess
calc_prob_coin_flip_data <- function(P_heads_guess, coin_flips)
    {
    # Empty list of probabilities
    probs_list = rep(NA, times=length(coin_flips))
    probs_list

    for (i in 1:length(coin_flips))
        {
        # Print an update
        #cat("\nAnalysing coin flip #", i, "/", length(coin_flips), sep="")

        # Get the current coin flip
        coin_flip = coin_flips[i]

        # If the coin flip is heads, give that datum
        # probability P_heads_guess.
        # If tails, give it (1-P_heads_guess)

        if (coin_flip == "H")
            {
            probs_list[i] = P_heads_guess
            } # End if heads

        if (coin_flip == "T")        
            {
            probs_list[i] = (1-P_heads_guess)
            } # End if tails
        } # End for-loop

    # Look at the resulting probabilities
    probs_list

    # We get the probability of all the data by multiplying
    # all the probabilities
    likelihood_of_data_given_P_heads_guess = prod(probs_list)

    # Return result
    return(likelihood_of_data_given_P_heads_guess)
    }

# Now, we can just use this function:
calc_prob_coin_flip_data(P_heads_guess=0.5, coin_flips=coin_flips)
calc_prob_coin_flip_data(P_heads_guess=0.6, coin_flips=coin_flips)
calc_prob_coin_flip_data(P_heads_guess=0.7, coin_flips=coin_flips)

# Look at that!  We did all of that work in a split-second.

# In fact, we can make another for-loop, and search for the ML
# value of P_heads by trying all of the values and plotting them.

# Sequence of 50 possible values of P_heads between 0 and 1
P_heads_values_to_try = seq(from=0, to=1, length.out=50)
likelihoods = rep(NA, times=length(P_heads_values_to_try))

for (i in 1:length(P_heads_values_to_try))
    {
    # Get the current guess at P_heads_guess
    P_heads_guess = P_heads_values_to_try[i]

    # Calculate likelihood of the coin flip data under
    # this value of P_heads
    likelihood = calc_prob_coin_flip_data(P_heads_guess=P_heads_guess, coin_flips=coin_flips)

    # Store the likelihood value
    likelihoods[i] = likelihood
    } # End for-loop

# Here are the resulting likelihoods:
likelihoods

# Let's try plotting the likelihoods to see if there's a peak
plot(x=P_heads_values_to_try, y=likelihoods)
lines(x=P_heads_values_to_try, y=likelihoods)

# Whoa! That's quite a peak!  You can see that the likelihoods
# vary over several orders of magnitude.
#
# Partially because of this extreme variation, we often use the 
# log-likelihood (natural log, here) instead of the raw
# likelihood.
#
# (Other reasons: machines have a minimum precision, log-likelihoods
#  can be added instead of multiplied, AIC is calculated from 
#  log-likelihood, etc.)
#
#
log_likelihoods = log(likelihoods, base=exp(1))

plot(x=P_heads_values_to_try, y=log_likelihoods)
lines(x=P_heads_values_to_try, y=log_likelihoods)

# Let's plot these together
par(mfrow=c(2,1))
plot(x=P_heads_values_to_try, y=likelihoods, main="Likelihood (L) of the data")
lines(x=P_heads_values_to_try, y=likelihoods)

plot(x=P_heads_values_to_try, y=log_likelihoods, main="Log-likelihood (LnL) of the data")
lines(x=P_heads_values_to_try, y=log_likelihoods)

# Maximum likelihood optimization
# 
# You can see that the maximum likelihood of the data occurs when 
# P_heads is somewhere around 0.6 or 0.7.  What is it 
# exactly?
#
# We could just keep trying more values until we find whatever
# precision we desire.  But, R has a function for
# maximum likelihood optimization!
# 
# It's called optim().  Optim() takes a function as an input.
# Fortunately, we've already written a function!
#
# Let's modify our function a bit to return the log-likelihood,
# and print the result:

# Function that calculates the probability of coin flip data
# given a value of P_heads_guess
calc_prob_coin_flip_data2 <- function(P_heads_guess, coin_flips)
    {
    # Empty list of probabilities
    probs_list = rep(NA, times=length(coin_flips))
    probs_list

    for (i in 1:length(coin_flips))
        {
        # Print an update
        #cat("\nAnalysing coin flip #", i, "/", length(coin_flips), sep="")

        # Get the current coin flip
        coin_flip = coin_flips[i]

        # If the coin flip is heads, give that datum
        # probability P_heads_guess.
        # If tails, give it (1-P_heads_guess)

        if (coin_flip == "H")
            {
            probs_list[i] = P_heads_guess
            } # End if heads

        if (coin_flip == "T")        
            {
            probs_list[i] = (1-P_heads_guess)
            } # End if tails
        } # End for-loop

    # Look at the resulting probabilities
    probs_list

    # We get the probability of all the data by multiplying
    # all the probabilities
    likelihood_of_data_given_P_heads_guess = prod(probs_list)

    # Get the log-likelihood
    LnL = log(likelihood_of_data_given_P_heads_guess)
    LnL

    # Error correction: if -Inf, reset to a low value
    if (is.finite(LnL) == FALSE)
        {
        LnL = -1000
        }

    # Print some output
    print_txt = paste("\nWhen P_heads=", P_heads_guess, ", LnL=", LnL, sep="")
    cat(print_txt)

    # Return result
    return(LnL)
    }

# Try the function out:
LnL = calc_prob_coin_flip_data2(P_heads_guess=0.1, coin_flips=coin_flips)
LnL = calc_prob_coin_flip_data2(P_heads_guess=0.2, coin_flips=coin_flips)
LnL = calc_prob_coin_flip_data2(P_heads_guess=0.3, coin_flips=coin_flips)

# Looks like it works!  Let's use optim() to search for he 
# best P_heads value:

# Set a starting value of P_heads
starting_value = 0.1

# Set the limits of the search
limit_bottom = 0
limit_top = 1

optim_result = optim(par=starting_value, fn=calc_prob_coin_flip_data2, coin_flips=coin_flips, method="L-BFGS-B", lower=limit_bottom, upper=limit_top, control=list(fnscale=-1))

# You can see the search print out as it proceeds.

# Let's see what ML search decided on:
optim_result

# Let's compare the LnL from ML search, with the binomial mean
optim_result$par

# Here's the formula:
P_heads_ML_estimate = numHeads / numTotal
P_heads_ML_estimate

# Wow!  Pretty good!  
#
# But -- why would anyone ever go through all the rigamarole, when they could 
# just calculate P_head directly?
#
# Well, only in simple cases do we have a formula for the maximum likelihood 
# estimation of the mean.  The optim() strategy works whether or not
# there is a simple formula.
#
# In real life science, ML optimization gets use A LOT, but most scientists
# don't learn it until graduate school, if then.
# 
# For a real-life example of ML analysis, try the tutorial for my biogeography
# R package, BioGeoBEARS:
# 
# http://phylo.wikidot.com/biogeobears#toc16
# 

# NOTE: BAYESIAN METHODSs
# By the way, having done this ML search, we are very close to being able 
# to do a Bayesian MCMC (Markov-Chain, Monte-Carlo) analysis.  However,
# we don't have time for this today.  Come talk to me this 
# summer if you are interested!

#######################################################
# CHAPTER 6 (BONUS): R PACKAGES AND MORE R
#
# Many good functions are found in base R, but there 
# are many, many more in other R packages
#######################################################

#############################################
# install a package (only have to do once)
#############################################
# Type this:
odd(13)

# What happened?
#

# Now do this:
library(gtools)							 # check if you have gtools	
# install.packages("gtools") # if needed
# gtools contains functions for odd/even (and many other things)

# after a package is installed, you have to library() it to use its
# functions during an R session
library(gtools)

# Now type this:
odd(13)

# For loops
# Here, we are using for loops, if statements, and the gtools function "odd"
# What does the code in this for loop do?
# 
for (i in 1:10)
    {
    if (odd(i) == TRUE)
        {
        print(paste(i, "is odd!", sep=" "))
        }
    else
        {
        print("Blah!")
        }
    }

#
paste("This", "is", "fun", sep=" ")

# print can be annoying, use cat
for (i in 1:10)
    {
    if (odd(i) == TRUE)
        {
        cat(paste(i, "is odd!", "\n",  sep=" "))
        }
    else
        {
        cat("Blah!\n" )
        }
    }

# How to make your own function
# (These can be sources from a source script, which is
#  a .R file either on your computer
get_square <- function(x)
    {
    output = x^2
    return(output)
    }

x = 4

#
(newval = get_square(x))

# Write to tab-delimited text file
fn = "grade_data.txt"
write.table(grade_data, file=fn, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

# read it back in
new_data = read.table(fn, header=TRUE, sep="\t", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)

# plots and stats

#
plot(new_data$grade1, new_data$grade2)

#
title("Hi Tom!")

#
plot(new_data$grade1, new_data$grade2, xlab="scores for grade #1", ylab="scores for grade #2")

#
lines(new_data$grade1, new_data$grade2)

#
pairs(new_data[, 2:4])

#
cor(new_data[, 2:4])

# The function "pairs_with_hist" is one I modified from 
# somewhere and put in genericR_v1, which lives at:
# https://gist.github.com/nmatzke/b21fcccb09b42959897d
# Get it with the source() command:
source("https://gist.githubusercontent.com/nmatzke/b21fcccb09b42959897d/raw/5f461b8eb1032fa7aad2e0ae5c13e9f9a1783cdb/_genericR_v1.R")
pairs_with_hist(new_data[, 2:4])

# CTRL-left or right to arrow through plots

# help on any function
?mean
?std

# once you are looking at the help page, go to the bottom & click index to see all the options for the package
# or just e.g. 
?gtools

# search for text in help
# (marginally useful)
??histogram

# I like to google: ' "r-help" something something '
# ...since someone has always asked my question on the r-help listserv

###############################################
# Basic crash course in APE 
# The R Package APE: Analysis of Phylogenetics
#   and Evolution
#
# Paradis's book on APE is linked from the 
# course website:
# http://ib.berkeley.edu/courses/ib200b/IB200B_SyllabusHandouts.shtml
# (for Feb. 3)
###############################################
# Install APE
library(ape) # check if you have ape
# install.packages("ape")
# (This should install some other needed packages also)

library(ape)

# This is what a Newick string looks like:
newick_str = "(((Humans, Chimps), Gorillas), Orangs);"
tr = read.tree(text=newick_str)
# If loading a tree file
# tr = read.tree(file=trfn)
plot(tr)

# What is the data class of "tr"?
#
class(tr)

# Is there any difference in the graphic produced by these two commands?
#
plot(tr)
plot.phylo(tr)

# What is the difference in the result of these two help commands?
#
?plot
?plot.phylo

# What are we adding to the tree and the plot of the tree, this time?
#
newick_str = "(((Humans:6.0, Chimps:6.0):1.0, Gorillas:7.0):1.0, Orangs:8.0):1.0;"
tr = read.tree(text=newick_str)
plot(tr)

# What are we adding to the tree and the plot of the tree, this time?
#
newick_str = "(((Humans:6.0, Chimps:6.0)LCA_humans_chimps:1.0, Gorillas:7.0)LCA_w_gorillas:1.0, Orangs:8.0)LCA_w_orangs:1.0;"
tr = read.tree(text=newick_str)
plot(tr, show.node.label=TRUE)

# More on Newick format, which, annoyingly, is sometimes inconsistent:
# http://en.wikipedia.org/wiki/Newick_format

# Have a look at how the tree is stored in R
tr

#
tr$tip.label

#
tr$edge

#
tr$edge.length

#
tr$node.label

# If you forget how to find these, you can use the "attributes" function
#
attributes(tr)



#######################################################
# Printing a tree table with BioGeoBEARS prt() function
#######################################################

trtable = prt(tr, printflag=FALSE)
trtable

# It can be useful to have a "universal identifier" for a node
# across different tree.
# Look at "tipnames" column:
trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE)

# If you have fossils, you might need to identify them automatically

tr_w_fossils_str = "(((Humans:6.0, Chimps:6.0)LCA_humans_chimps:1.0, Gorillas:7.0)LCA_w_gorillas:1.0, Orangs:4.0)LCA_w_orangs:1.0;"
tr3 = read.tree(text=tr_w_fossils_str)
plot(tr3)
trtable = prt(tr3, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=0.5)
trtable



# Now plot the tree in different ways:
# (CTRL-right or CTRL-left to flip between the trees in the graphics window)
plot(tr, type="phylogram")

#
plot(tr, type="phylogram", direction="rightwards")
plot(tr, type="phylogram", direction="leftwards")
plot(tr, type="phylogram", direction="upwards")
plot(tr, type="phylogram", direction="downwards")

#
plot(tr, type="cladogram")
plot(tr, type="fan")
plot(tr, type="unrooted")
plot(tr, type="radial")

#
plot(tr, type="unrooted", edge.width=5)
plot(tr, type="unrooted", edge.width=5, edge.color="blue")
plot(tr, type="unrooted", edge.width=5, edge.color="blue", lab4ut="horizontal")
plot(tr, type="unrooted", edge.width=5, edge.color="blue", lab4ut="axial")

# In R GUI, you can save any displayed tree to PDF, or do a screen capture etc.
# you can also save a tree to PDF as follows:
pdffn = "homstree.pdf"
pdf(file=pdffn)
plot(tr, type="unrooted", edge.width=5, edge.color="blue", lab4ut="axial")
dev.off()

# In Macs (and maybe PCs), this will open the PDF from R:
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

# How to save the tree as text files
#
newick_fn = "homstree.newick"
write.tree(tr, file=newick_fn)

#
nexus_fn = "homstree.nexus"
write.nexus(tr, file=nexus_fn)


# To conclude the lab, I wanted to find, download, and display
# a "tree of life".
#
# To do this, I went to the TreeBase search page:
# http://www.treebase.org/treebase-web/search/studySearch.html
# (NOTE 2019: TreeBase now gives a 404 error!)d
#
# ...and searched on studies with the title "tree of life"
#
# Annoyingly, the fairly famous tree from:
#
# Ciccarelli F.D. et al. (2006). "Toward automatic reconstruction of 
# a highly resolved tree of life." Science, 311:1283-1287.
# http://www.sciencemag.org/content/311/5765/1283.abstract
#
# ...was not online, as far as I could tell.  And a lot of these are the "turtle trees of life", etc.
# Lame.  But this one was a tree covering the root of known
# cellular life.
#
# Caetano-anollés G. et al. (2002). "Evolved RNA secondary structure
# and the rooting of the universal tree of life." Journal of
# Molecular Evolution.
#
# I got tree S796 for this study, then click over to the "Trees" tab to get the
# tree...
#
# http://www.phylowidget.org/full/?tree=%27http://www.treebase.org/treebase-web/tree_for_phylowidget/TB2:Tr3931%27
# (2019: link also doesn't work!)
#
# Or, download the tree from our website, here:
# http://ib.berkeley.edu/courses/ib200b/labs/Caetano-anolles_2002_JME_ToL.newick
#
# Or, click on "Files", at the bottom-right of 
# http://phylo.wikidot.com/introduction-to-r-pcms
#
# I.e.: 

# load the tree and play with it:
newick_fn = "Caetano-anolles_2002_JME_ToL.newick"
tree_of_life = read.tree(newick_fn)
plot(tree_of_life, type="cladogram")
plot(tree_of_life, type="phylogram")
plot(tree_of_life, type="unrooted", lab4ut="axial")

# aw, no branch lengths in TreeBase! Topology only! Lame!




















#######################################################
# LAST ADDITION - a few helpful BioGeoBEARS functions
# for working with trees:
#######################################################
# if you have this installed:
library(BioGeoBEARS)

# To install BioGeoBEARS, see:
# https://github.com/nmatzke/BioGeoBEARS

# "moref" -- the equivalent of the bash shell "more", in R:
nexus_fn = "homstree.nexus"
moref(nexus_fn)





# "prt" -- prints a tree to a table data.frame, in R:

newick_fn = "homstree.newick"
tr = read.tree(file=newick_fn)
prt(tr)


# Alphabetical, concatenated tipnames as universal unique node labels:
prt(tr, get_tipnames=TRUE)

# Check how fossils are labeled:
prt(tr, fossils_older_than=0.01, get_tipnames=TRUE)









# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size = 4
#areas = getareas_from_tipranges_object(tipranges)
areas = c("K", "O", "M", "H")

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based))
{    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list