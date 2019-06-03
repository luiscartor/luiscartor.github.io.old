R Mini Tutorial: R Environment, reading/writting, and basic executions
================
luiscartor
May 29, 2019

2.1 The R environment
---------------------

You should have installed and be familiar with R and RStudio by now.

In R, we talk about the **workspace** as your current R working environment and includes any user-defined objects (vectors, matrices, data frames, lists, functions). At the end of an R session, the user can save an image of the current workspace that is automatically reloaded the next time R is started.

### 2.1.1 Working directory

By setting the working directory of R, we can directly use relative path to read/write data within this folder.

For example: an absolute path looks like this:
`D:/SRE2019/Rtutorial/part2/exercise1/occurrence.csv`
if we set the working directory as `D:/SRE2019/Rtutorial/part2`, then the relative path will be:
`exercise1/occurrence.csv`

``` r
setwd("D:/SRE2019/Rtutorial/part2") # note: the folder should exist, before you set it up.
```

In RStudio, you can set up the working directory in Session&gt;Set Working Directory.

You can check what your actual working directory is by using the `getwd` command:

``` r
getwd()
```

(Note that if you are using a Microsoft system the file naming convention is different from what we use here. If you want to use a backslash it needs to be escaped, i.e. use two backslashes together)

### 2.1.2 Installing and loading packages

We often need additional functionality beyond those offered by the core R library. In order to install an extension package, you should use install.packages function at the prompt and follow the instructions:

``` r
install.packages("raster") # package for raster analysis
```

It is also possible to use RStudio to install packages. Just go to Tools &gt; Install Packages. Or open the "Packages" tab on the right-bottom window and click Install.

Now we need to load the package into our session:

``` r
library("raster")
```

For the next lesson on spatial analysis, we will need a bunch of packages. We can install several packages at the same time by creating a vector with each package name as an element.

``` r
packages_needed <- c("dismo", # a collection of ENM/SDM tools
                     "rgeos","rgdal","sp", # spatial data analysis
                     "ENMeval", # a few new tools in ENM/SDM
                     "wallace",   # interface for Maxent and ENMeval
                     "utils", # for zip & unzip files
                     "jsonlite" # necessary for download data from GBIF
                     )
pk_to_install <- packages_needed [!( packages_needed %in% rownames(installed.packages())  )]
if(length(pk_to_install)>0 ){
  install.packages(pk_to_install,repos="http://cran.r-project.org")
}
library("raster")
library("dismo")
library("rgeos")
library("rgdal")
library("sp")
library("ENMeval")
```

Note that we used an `if` statement in the last piece of code. This is telling R to install the packages, only if they hadn't been installed before. We will study the "if" statement in a few minutes.

<br>

> ### Excercise
>
> Install and load the graphics package *ggplot2*.

<br>

### 2.1.3 Getting help

R includes a comprehensive built-in help system with extensive documentation. You can access the help from the command prompt. Try the following commands:

``` r
help.start()   # general help
```

    ## starting httpd help server ... done

    ## If the browser launched by 'xdg-open' is already running, it is
    ##     *not* restarted, and you must switch to its window.
    ## Otherwise, be patient ...

``` r
help(lm)      # help about function lm (linear models)
?anova           # same thing as above
example(lm)   # show an example of function anova
```

    ## 
    ## lm> require(graphics)
    ## 
    ## lm> ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
    ## lm> ## Page 9: Plant Weight Data.
    ## lm> ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
    ## 
    ## lm> trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
    ## 
    ## lm> group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
    ## 
    ## lm> weight <- c(ctl, trt)
    ## 
    ## lm> lm.D9 <- lm(weight ~ group)
    ## 
    ## lm> lm.D90 <- lm(weight ~ group - 1) # omitting intercept
    ## 
    ## lm> ## No test: 
    ## lm> ##D anova(lm.D9)
    ## lm> ##D summary(lm.D90)
    ## lm> ## End(No test)
    ## lm> opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
    ## 
    ## lm> plot(lm.D9, las = 1)      # Residuals, Fitted, ...

![](RMT_part2_files/figure-markdown_github/unnamed-chunk-5-1.png)

    ## 
    ## lm> par(opar)
    ## 
    ## lm> ## Don't show: 
    ## lm> ## model frame :
    ## lm> stopifnot(identical(lm(weight ~ group, method = "model.frame"),
    ## lm+                     model.frame(lm.D9)))
    ## 
    ## lm> ## End(Don't show)
    ## lm> ### less simple examples in "See Also" above
    ## lm> 
    ## lm> 
    ## lm>

If you are not sure about the name of the function you are looking for, you can perform a fuzzy search with the apropos function:

``` r
apropos("anova")
```

    ## [1] "anova"            "manova"           "power.anova.test"
    ## [4] "stat.anova"       "summary.manova"

``` r
help(anova)
```

<br>

> ### Excercise
>
> Explore the documentation and examples of the *ggplot2* library.

<br>

2.2 Reading and writing
-----------------------

In R, we can import/read data stored outside the R environment. We can also write data into files which will be stored and accessed by the operating system. R offers options to import many file types formats like txt, csv, excel, etc. Additional packages will allow us to read and write other types of data, such as spatial data.

### 2.2.1 Reading files

Let's import and read a CSV file. First open your spread sheet software and create a simple table. Give names to the columns and give characters and numerical values. For example:

| heron      | colony | population |
|------------|--------|------------|
| garzetta   | 1      | 55         |
| garzetta   | 2      | 78         |
| nycticorax | 1      | 25         |
| nycticorax | 3      | 60         |
| alba       | 1      | 36         |
| cinerea    | 4      | 80         |

Name it herons.csv and save it in your working directory.

Now we read the file in R:

``` r
herons <- read.csv(file="herons.csv",head=TRUE,sep=",") # first row contains variable names, comma is separator 
herons
```

    ##        heron colony population
    ## 1   garzetta      1         55
    ## 2   garzetta      2         78
    ## 3 nycticorax      1         25
    ## 4 nycticorax      3         60
    ## 5       alba      1         36
    ## 6    cinerea      4         80

We can summarize the new table:

``` r
summary(herons)
```

    ##         heron       colony       population   
    ##  alba      :1   Min.   :1.00   Min.   :25.00  
    ##  cinerea   :1   1st Qu.:1.00   1st Qu.:40.75  
    ##  garzetta  :2   Median :1.50   Median :57.50  
    ##  nycticorax:2   Mean   :2.00   Mean   :55.67  
    ##                 3rd Qu.:2.75   3rd Qu.:73.50  
    ##                 Max.   :4.00   Max.   :80.00

You can explore each column using the subsetting tools that you already learn in Part 1:

``` r
herons$population
```

    ## [1] 55 78 25 60 36 80

``` r
herons["alba",]
```

    ##    heron colony population
    ## NA  <NA>     NA         NA

### 2.2.2 Writing files

Before learning how to write a file, let's edit the previous talbe. Let's add an extra row:

``` r
herons_edit <- rbind(herons, c("alba", 5, 12))
```

Now we write our edited table:

``` r
write.csv(herons_edit,"herons_edit.csv",row.names = FALSE)
```

Check that the table has been written to disk and import it again to R:

``` r
new_table <- read.csv(file="herons_edit.csv",head=TRUE,sep=",")
new_table
```

    ##        heron colony population
    ## 1   garzetta      1         55
    ## 2   garzetta      2         78
    ## 3 nycticorax      1         25
    ## 4 nycticorax      3         60
    ## 5       alba      1         36
    ## 6    cinerea      4         80
    ## 7       alba      5         12

<br>

> ### Excercise
>
> Use the function `read.table()` to import the herons table. Are there any differences between `read.table()` and `read.csv()`? Use help() to assist you in using this command.

<br>

2.3 Basic executions
--------------------

### 2.3.1 Conditional execution

### 2.3.2 Loops

2.4 Statistical models
----------------------

### 2.4.1 Linear models

### 2.4.2 ANOVAs

2.5 Plotting in R
-----------------

### 2.5.1 The plot() function
