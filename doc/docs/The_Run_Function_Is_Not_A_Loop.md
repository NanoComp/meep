---
The Run Function Is Not A Loop
---

In Meep, there are functions `run-until` and similar that run the simulation, and take arguments allowing custom actions to be performed on every time step, or on some subset of the time steps.  Many users misunderstand this, however, and make the same mistake: they think the `run` function is a "looping" construct of some kind, and that you can just put any code you want into it and it will get executed for every time step.  This mistake and how to correct it are described in this article.

Hello World
-----------

Let's consider a "Hello World" example.  Suppose we start with a control file that runs for 200 time units and outputs <math>E_z</math> on each time step:

```scm
 (run-until 200 output-efield-z)
```

and now we want to modify it to also print "Hello World!" for every time step, as it is running.

The Wrong Way
-------------

Many users will naively write:

```scm
 (run-until 200 output-efield-z (print "Hello World!\n"))     
```

**This is wrong.**  It will output "Hello World!" *once*, then give an error.  What is going on? The problem is that you are thinking of `run-until` in the wrong way, as if it were a loop:

```scm
  for time < 200 do
      output-efield-z
      (print "Hello World!\n")
```

This is **not** what is happening. Instead, `run-until` is just a *function* that runs the simulation, and its arguments should be *functions* that are called for each time step.  That is, it is really doing something like:

```html
  evaluate the arguments: 200: a number
           output-efield-z: a function
           (print "Hello World!\n"): prints output and returns #<unspecified>;
  call run-until:
     time-step until t=200
     at each time step, call the arguments:
        call (output-efield-z)
        call (#<unspecified>)
```

Two things went wrong.  First, the arguments are evaluated **before** calling the function, which means that the `print` statement is executed before `run-until` even starts!  Second, `run-until` then tries to call the **result** of `(print ...)` as if it were a function, which causes an error because `(print ...)` does not return a function. The `print` returns a special Scheme code `#<unspecified>` that means it doesn't really return anything at all.

The Right Way
-------------

What we should have passed to `run-until`, instead of the *result* of calling `(print ...)`, is a *function* that calls `(print ...)`.  There are two ways to do it.

First, we could explicitly define a function, call it `my-hello`, that does what we want:

```scm
 (define (my-hello) (print "Hello World!\n"))
 (run-until 200 output-efield-z my-hello)
```

Notice two things. First, `my-hello` is a function of no arguments, which means that it is just called at every time step. Another, more complicated, possibility is described in the [User Interface](Scheme_User_Interface.md#writing-your-own-step-functions). Second, when we call `run-until`, we just pass the *name* of the function `my-hello`, and not the *result* of calling the function `(my-hello)`.

A second possibility is that we could use Scheme's `lambda` construct to define our function in-line. The `lambda` syntax in Scheme allows you to define *anonymous* functions without assigning them a name via `define`, and to stick the function definition right into another expression. It works like this:

```scm
 (run-until 200 output-efield-z (lambda () (print "Hello World!\n")))
```

Here, the `(lambda () ...)` defines a function of no arguments `()` that, when called, executes the `...` statements.