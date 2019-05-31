R Mini Tutorial: Getting started and data types
================
luiscartor
May 29, 2019

1. Getting started
------------------

In Unix or Mac, execute R on the command line, or click on the R icon on the desktop for Windows.

To finish an R session, type `q()` at the command prompt. You will be asked to save the workspace when finishing the session. It is a good practice to NOT to save the session, but if you say yes, the objects on your workspace will be saved as a .RData file. If saved, these workspace will be loaded when you open a new session.

1.1 Basic calculations
----------------------

Use the console to execute R commands, such as basic calculations. After hitting enter, the console will show the result.

``` r
  3+2
```

    ## [1] 5

But the result of the calculation can be also stored to a variable:

``` r
  a <- 3+2
```

Variable names in R must begin with a letter, followed by alphanumeric characters. Periods or underscores are often used to break up long variables (e.g. `long_variable`). Use variable names that help you identify the variable, but without being very long (e.g. herons\_popmean). Avoid `c`, `q`, `t`, `C`, `D`, `F`, `I`, `T`, and things like `mean`, `var`, `sd`, as they are built-in R functions and they will cause errors.

Now '5' is stored in the variable a. To print the value of a variable, just type the variable name:

``` r
  a
```

    ## [1] 5

Some arithmetic operators:

| Operator  | Description    |
|-----------|----------------|
| +         | addition       |
| -         | subtraction    |
| \*        | multiplication |
| /         | division       |
| ^ or \*\* | exponentiation |

Some comparison and logical operators:

| Operator | Description              |
|----------|--------------------------|
| &gt;     | greater than             |
| &gt;=    | greater than or equal to |
| ==       | exactly equal to         |
| !=       | not equal to             |
| !        | logical NOT              |
| &        | logical AND              |

Now we can start creating variables and making arithmetic and logical operations.

``` r
x <- 5
y <- 2
z1 <- x*y
z2 <- x/y
z3 <- x^y
z1; z2; z3
```

    ## [1] 10

    ## [1] 2.5

    ## [1] 25

Note that we can use a semicolon to put two or more commands in one line. In the last line of code, we use the semicolon to print the values of the three z variables.

We can combine several operations:

``` r
A <- 5; B <- 3
C <- (A+B)/(A+5*B) 
C
```

    ## [1] 0.4

And use the logical operators:

``` r
A == B
```

    ## [1] FALSE

``` r
A > B
```

    ## [1] TRUE

``` r
C != B
```

    ## [1] TRUE

We can use a wide rage of built-in mathematical functions in R

| R command              | function                          |
|------------------------|-----------------------------------|
| abs()                  | absolute value, $                 |
| cos(), sin(), tan()    | cosine, sine, tangent             |
| acos(), asin(), atan() | arc-cosine, arc-sine, arc-tangent |
| exp(x)                 | exponential function              |
| log()                  | natural (base-e) logarithm        |
| log10()                | base-10 logarithm                 |
| sqrt()                 | square root, √x                   |

For example, try this calculations:

``` r
A <- 5; B <- 3
C <- (A+2*sqrt(B))/(A+5*sqrt(B)) 
D <- log(A+B)/cos(C)
C;D
```

    ## [1] 0.6196152

    ## [1] 2.554277

As usual, parentheses alter the order of the operations. Try:

``` r
C <- A+2*sqrt(B)/A+5*sqrt(B)
C
```

    ## [1] 14.35307

and

``` r
C <- A + 2*(sqrt(B)/A + 5) *sqrt(B)
C
```

    ## [1] 23.52051

<br>

> ### Excercise: Solve the following equations
>
> 1.  
>     1 + 0.3 + 0.7<sup>2</sup>
<<<<<<< HEAD
> 2.  For a = 0.5 and b = 1:
>     *s**i**n*(*a*)*c**o**s*(*b*)−*c**o**s*(*a*)
> 3.  
=======
<<<<<<< HEAD
>      &lt;&lt;&lt;&lt;&lt;&lt;&lt; HEAD
> 2.  For a = 0.5 and b = 1: 
>     *s**i**n*(*a*)*c**o**s*(*b*)−*c**o**s*(*a*)
>     ===========================================
>
> 3.  For a = 0.5 and b = 1:
>     $$\\frac{sin(a)cos(b)}{cos(a)}$$
>      &gt;&gt;&gt;&gt;&gt;&gt; 7ee0d78e0545a65101691046db9156a8f4cf96ac
> 4.  
=======
> 2.  For a = 0.5 and b = 1:
>     *s**i**n*(*a*)*c**o**s*(*b*)−*c**o**s*(*a*)
> 3.  
>>>>>>> 71ba72a3ddbdca1342d4a3ca81e9f0a690378873
>>>>>>> 25774a4b398458984a2042ad37f3ec89908c1986
>     *e*<sup>0.5<sup>2</sup>/4</sup>
> 4.  For a = 10 and b = 2:
>     log(*a*)−3 \* log 10(*b*)

<br>

1.2 Basic Data Types
--------------------

### 1.2.1 Vectors

The simplest data type in R is the vector, which can be: numerical, character, logical, factor, or list based.

``` r
  a <- c(1,2,3,7.3,-2,-0.1) # numeric vector
  b <- c("one","two","three","a","b","c") # character vector
  c <- c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE) #logical vector
  a;b;c
```

    ## [1]  1.0  2.0  3.0  7.3 -2.0 -0.1

    ## [1] "one"   "two"   "three" "a"     "b"     "c"

    ## [1]  TRUE  TRUE  TRUE FALSE  TRUE FALSE

We can access certain elements of the vector by using subscripts:

``` r
a[5]
```

    ## [1] -2

``` r
b[3]
```

    ## [1] "three"

``` r
a[c(1,3,5)]
```

    ## [1]  1  3 -2

Arithmetic expressions can be used with numerical vectors:

``` r
  a <- a+1; a
```

    ## [1]  2.0  3.0  4.0  8.3 -1.0  0.9

``` r
  a <- sqrt(a); a
```

    ## Warning in sqrt(a): NaNs produced

    ## [1] 1.4142136 1.7320508 2.0000000 2.8809721       NaN 0.9486833

Also expressions involving multiple vectors:

``` r
  a <- c(1,2,3,7.3,-2,-0.1)
  b <- c(1,0,1,0,2,2.5)
  c <- 2*a + b + sqrt(c) + 8
  c
```

    ## [1] 12.0 13.0 16.0 22.6  7.0 10.3

What happens when the vectors have different lengths? The elements of the shorter vector are re-utilized for the elements of the long vector:

``` r
  a <- c(10,20,30,40,50,60)
  b <- c(1,2,3)
  a + b
```

    ## [1] 11 22 33 41 52 63

Some basic functions can be used to create vectors. One of this functions is `seq`, with syntax: `a <- seq(from, to, by)`:

``` r
  a <- seq(1,9,2)
  a
```

    ## [1] 1 3 5 7 9

Another useful function to create vectors is `rep` (`b <- rep(values,length)`):

``` r
  b <- rep(3,5)
  b
```

    ## [1] 3 3 3 3 3

### 1.2.2 Matrices

Now we are going to take a look how to create and handle matrices in R. All columns in a matrix must have the same data mode (numeric, character, logical, etc.) and the same length. A matrix is created using the following syntax: `m <- matrix(vector, nrow=r, ncol=c, byrow=FALSE, dimnames=list(rnames, cnames))`. See the example:

``` r
  m <- matrix(1:20, nrow=5, ncol=4)
  m
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    1    6   11   16
    ## [2,]    2    7   12   17
    ## [3,]    3    8   13   18
    ## [4,]    4    9   14   19
    ## [5,]    5   10   15   20

Elements of the matrix can be accessed using subscripts:

``` r
  m[,3] # 3th column of matrix
```

    ## [1] 11 12 13 14 15

``` r
  m[4,] # 4rd row of matrix
```

    ## [1]  4  9 14 19

``` r
  m[2:4,1:3] # rows 2,3,4 of columns 1,2,3 
```

    ##      [,1] [,2] [,3]
    ## [1,]    2    7   12
    ## [2,]    3    8   13
    ## [3,]    4    9   14

A matrix can be also created by providing a vector and the row and column names:

``` r
  v <- c(1,12,24,36)
  rnames <- c("R1", "R2") # character vector with row names
  cnames <- c("C1", "C2") # character vector with column names
  m <- matrix(v, nrow=2, ncol=2, byrow=TRUE,
  dimnames=list(rnames, cnames)) 
  m
```

    ##    C1 C2
    ## R1  1 12
    ## R2 24 36

An array is another type of matrix, but it can have more than two dimensions:

``` r
v1 <- c(5,9,3)
v2 <- c(10,11,12,13,14,15)

arr <- array(c(v1,v2),dim = c(3,3,2))
arr
```

    ## , , 1
    ## 
    ##      [,1] [,2] [,3]
    ## [1,]    5   10   13
    ## [2,]    9   11   14
    ## [3,]    3   12   15
    ## 
    ## , , 2
    ## 
    ##      [,1] [,2] [,3]
    ## [1,]    5   10   13
    ## [2,]    9   11   14
    ## [3,]    3   12   15

There are two very useful functions to deal with matrices, and that help us combine vectors of same size intro matrices, or combine matrices with vectors: `cbind` and `rbind`:

``` r
  m1 <- cbind(1:3,4:6,5:7) 
  m1
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    4    5
    ## [2,]    2    5    6
    ## [3,]    3    6    7

``` r
  m2 <- rbind(1:3,4:6)
  m2
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    2    3
    ## [2,]    4    5    6

<br>

> ### Excercises
>
> 1.  Try to combine matrices m1 and m2 into one single matrix. What are the dimensions of the new matrix?
> 2.  Subset the two first rows of m1 and combine the resulting matrix with m2. What are the dimensions of the new matrix?

<br>

### 1.2.2 Data frames

A data frame is similar to a matrix, but different columns can have different data modes (numeric, character, factor, etc.).

``` r
  a <- c(1,2,3,7.3,-2,-0.1) # numeric vector
  b <- c("one","two","three","a","b","c") # character vector
  c <- c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE) #logical vector
  df <- data.frame(a,b,c)
  names(df) <- c("numbers","strings","passed")
  df
```

    ##   numbers strings passed
    ## 1     1.0     one   TRUE
    ## 2     2.0     two   TRUE
    ## 3     3.0   three   TRUE
    ## 4     7.3       a  FALSE
    ## 5    -2.0       b   TRUE
    ## 6    -0.1       c  FALSE

We can use different methods to identify elements of the data frame:

``` r
df[1:2] # columns 2 and 3 of data frame
```

    ##   numbers strings
    ## 1     1.0     one
    ## 2     2.0     two
    ## 3     3.0   three
    ## 4     7.3       a
    ## 5    -2.0       b
    ## 6    -0.1       c

``` r
df[c("numbers","passed")] # columns numbers and passed 
```

    ##   numbers passed
    ## 1     1.0   TRUE
    ## 2     2.0   TRUE
    ## 3     3.0   TRUE
    ## 4     7.3  FALSE
    ## 5    -2.0   TRUE
    ## 6    -0.1  FALSE

``` r
df$strings # variable strings in the data frame 
```

    ## [1] one   two   three a     b     c    
    ## Levels: a b c one three two

### 1.2.3 Lists

A list can contain anything in each element: vectors, matrices, other lists and other objects:

``` r
L <- list(A=a,B=b,C=c)
L
```

    ## $A
    ## [1]  1.0  2.0  3.0  7.3 -2.0 -0.1
    ## 
    ## $B
    ## [1] "one"   "two"   "three" "a"     "b"     "c"    
    ## 
    ## $C
    ## [1]  TRUE  TRUE  TRUE FALSE  TRUE FALSE

To access an element of the list, we need to use `[[ ]]` by number, or `$` to access an element by name:

``` r
L[[1]] # Access elements A of the list by number
```

    ## [1]  1.0  2.0  3.0  7.3 -2.0 -0.1

``` r
L[c("B","C")] # Access elements B and C of the list (access by element name)
```

    ## $B
    ## [1] "one"   "two"   "three" "a"     "b"     "c"    
    ## 
    ## $C
    ## [1]  TRUE  TRUE  TRUE FALSE  TRUE FALSE

### 1.2.4 Factors

A factor is a variable that can take a finite number of distinct levels. So it is a **nominal** variable. We can build a factor by applying a `factor` function to a vector of any class:

``` r
v1 <- rep(c(1,2,3),each=3)
factor(v1)
```

    ## [1] 1 1 1 2 2 2 3 3 3
    ## Levels: 1 2 3

| Useful functions             | Description                      |
|------------------------------|----------------------------------|
| length(object)               | number of elements or components |
| str(object)                  | structure of an object           |
| class(object)                | class or type of an object       |
| names(object)                | names                            |
| c(object,object,...)         | combine objects into a vector    |
| cbind(object, object, ...)   | combine objects as columns       |
| rbind(object, object, ...)   | combine objects as rows          |
| object                       | prints the object                |
| ls()                         | list current objects             |
| rm(object)                   | delete an object                 |
| newobject &lt;- edit(object) | edit copy and save as newobject  |
| fix(object)                  | edit in place                    |

<br>

> ### Challenges
>
> 1.  Try to access only the string "three" from the list L created previously.
> 2.  Check the class type of each element of list L.
