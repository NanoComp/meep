---
# Guile and Scheme Information
---

There are many places you can go to find out more regarding Guile and the Scheme programming language. We list a few of them here.

[TOC]

Scheme
-------

Scheme is a simplified derivative of [Lisp](https://en.wikipedia.org/wiki/Lisp), and is a small and beautiful dynamically typed, [lexically scoped](https://en.wikipedia.org/wiki/Lexical_variable_scoping), [functional](https://en.wikipedia.org/wiki/Functional_programming_language) language.

-   A [history and introduction to Scheme](https://en.wikipedia.org/wiki/Scheme_programming_language)
-   [R5RS](http://www.swiss.ai.mit.edu/ftpdir/scheme-reports/r5rs-html/r5rs_toc.html) is the official Scheme language definition and reference.
-   A classic [introduction](ftp://ftp.cs.indiana.edu/pub/scheme-repository/doc/pubs/intro.txt) to Scheme by Ken Dickey.
-   [Structure and Interpretation of Computer Programs](http://mitpress.mit.edu/sicp/sicp.html) by Abelson, Sussman, and Sussman (full text online).
-   [Introduction to Scheme and its Implementation](ftp://ftp.cs.utexas.edu/pub/garbage/cs345/schintro-v14/schintro_toc.html) (the complete book on-line) by Prof. Paul R. Wilson ([Univ. of Texas](http://www.cs.utexas.edu/)).
-   [Teach Yourself Scheme](http://ds26gte.github.io/tyscheme/index.html) is a nice tutorial-style introduction to Scheme programming.
-   The [MIT Scheme Home Page](http://www.swiss.ai.mit.edu/projects/scheme/index.html) (where do you think Scheme was invented?)
    -   also check out the MIT [Scheme Underground](http://www.ai.mit.edu/projects/su/su.html)
-   There is the [comp.lang.scheme](news:comp.lang.scheme) newsgroup, and its [FAQ](http://www.faqs.org/faqs/by-newsgroup/comp/comp.lang.scheme.html).
-   The [Internet Scheme Repository](http://www.cs.indiana.edu/scheme-repository/) has a lot of code and documentation.
-   [schemers.org](http://www.schemers.org/) is another Scheme site and collection of resources.

Guile
------

Guile is an open-source implementation of Scheme, designed to be plugged in to other programs as a scripting language.

-   The [homepage](http://www.gnu.org/software/guile/) for the GNU Guile project.
-   See parts IV and V of the [Guile Reference Manual](http://www.gnu.org/software/guile/manual/html_node/index.html) for additional Scheme functions and types defined within the Guile environment.

How to Write a Loop in Scheme
-----------------------------

The most frequently asked question seems to be: **how do I write a loop in Scheme?** We give a few answers to that here, supposing that we want to vary a parameter *x* from *a* to *b* in steps of *dx*, and do something for each value of *x*.

The classic way, in Scheme, is to write a [tail-recursive](https://en.wikipedia.org/wiki/Tail_call) function:

`(define (doit x x-max dx)`
`   (if (<= x x-max)`
`      (begin`
`         `*`...perform` `loop` `body` `with` `x...`*
`         (doit (+ x dx) x-max dx))))`
`(doit a b dx) ; execute loop from a to b in steps of dx`

There is also a [do-loop construct](http://www.swiss.ai.mit.edu/ftpdir/scheme-reports/r5rs-html/r5rs_6.html#SEC36) in Scheme that you can use

`(do ((x a (+ x dx))) ((> x b)) `*`...perform` `loop` `body` `with` `x...`*`)`

If you have a list of values of *x* that you want to loop over, then you can use `map`:

`(map (lambda (x) `*`...do` `stuff` `with` `x...`*`) `*`list-of-x-values`*`)`

How to Read In Values from a Text File in Scheme
------------------------------------------------

A simple command to read a text file and store its values within a variable in Scheme is `read`. As an example, suppose a file *foo.dat* contains the following text, including parentheses:

`(1 3 12.2 14.5 16 18)`

In Scheme, we would then use

`(define port (open-input-file "foo.dat"))`
`(define foo (read port))`
`(close-input-port port)`

The variable *foo* would then be a list of numbers '(1 3 12.2 14.5 16 18).

Libctl Tricks Specific to [Meep](index.md) and [MPB](https://mpb.readthedocs.io)
--------------------------------------------------------------------------------

[libctl](https://libctl.readthedocs.io) has a couple of built-in functions `arith-sequence` and `interpolate` (see the [User Reference](https://libctl.readthedocs.io/en/latest/Libctl_User_Reference/)) to construct lists of a regular sequence of values, which you can use in conjunction with `map` as above:

`(map (lambda (x) `*`...do` `stuff` `with` `x...`*`) (arith-sequence x-min dx num-x))`

or

`(map (lambda (x) `*`...do` `stuff` `with` `x...`*`) (interpolate num-x (list a b)))`

Finally, if you have an entire libctl input file `myfile.ctl` that you want to loop, varying over some parameter *x*, you can do so by writing a loop on the Unix command-line. Using the [bash](https://en.wikipedia.org/wiki/bash) shell, you could do:

``for x in `seq a dx b`; do meep x=$x myfile.ctl; done``