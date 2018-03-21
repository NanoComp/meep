---
# Field Functions
---

As described in the [User Interface](Python_User_Interface.md), Meep provides several routines to integrate, analyze, and output arbitrary user-specified functions of the field components. See the functions whose names end with `_field_function`. This facility, while powerful, requires a bit more programming than most Meep usage, and is best illustrated by a few examples.

Every field-function that can be passed to these routines is of the form *f*(**r**,components...), where **r** is a position vector and "components..." are zero or more field components that the function depends on. The set of desired components is user-specified. As an example, suppose we are interested in the arbitrary function:

$$f(\mathbf{r}, E_x, H_z, \varepsilon) = x |\mathbf{r}| + E_x - \varepsilon H_z$$

We would define this function by:

*Python*
```py
def f(r, ex, hz, eps):
    return (r.x * r.norm() + ex) - (eps * hz)
```

*Scheme*
```scm
(define (f r ex hz eps)
   (- (+ (* (vector3-x r) (vector3-norm r)) ex) (* eps hz)))
```

Note that the `r` argument is a `Vector3` (Python) or `vector3` (Scheme), and can be manipulated by the functions defined in the [Libctl User Reference](https://libctl.readthedocs.io/en/latest/Libctl_User_Reference/).

Now, suppose we want to compute the integral of this function, over the whole computational cell. We can do this by calling the function `integrate_field_function` (Python) or `integrate-field-function` (Scheme), as follows:

*Python*
```py
print("The integral of our weird function is: {}"
	   .format(meep.Simulation.integrate_field_function([meep.Ex, meep.Hz, meep.Dielectric], f)))
```

*Scheme*
```scm
(print "The integral of our weird function is: "
       (integrate-field-function (list Ex Hz Dielectric) f) "\n")
```

Note that the first argument to `integrate_field_function` (Python) or `integrate-field-function` (Scheme) is a list (a standard Python/Scheme type) of `component` constants, specifying in order the list of field components our function `f` expects to be passed. Meep will then call `f` for every point in the computational cell in parallel on a parallel machine, and return the integral approximated by a [trapezoidal rule](https://en.wikipedia.org/wiki/trapezoidal_rule).

You can also specify an optional third argument to `integrate_field_function` or `integrate-field-function`, specifying an integration volume in case you don't want the integral over the whole computational cell. For example, the following code computes the integral of `f` along a line from (-1,0,0) to (1,0,0):

*Python*
```py
print("The integral of our weird function from (-1,0,0) to (1,0,0) is: {}"
	   .format(meep.Simulation.integrate_field_function([meep.Ex, meep.Hz, meep.Dielectric], f, meep.Volume(size=meep.Vector3(1), center=meep.Vector3()))))
```

*Scheme*
```scm
(print "The integral of our weird function from (-1,0,0) to (1,0,0) is: "
       (integrate-field-function (list Ex Hz Dielectric) f (volume (size 1 0 0) (center 0 0 0))) "\n")
```

Instead of computing the integral, Meep also provides a function to compute the maximum absolute value of our given function:

*Python*
```py
print("The maximum absolute value of our weird function from (-1,0,0) to (1,0,0) is: {}"
	   .format(meep.Simulation.max_abs_field_function([meep.Ex, meep.Hz, meep.Dielectric], f, meep.Volume(size=meep.Vector3(1), center=meep.Vector3()))))
```

*Scheme*
```scm
(print "The maximum absolute value of our weird function from (-1,0,0) to (1,0,0) is: "
       (max-abs-field-function (list Ex Hz Dielectric) f (volume (size 1 0 0) (center 0 0 0))) "\n")
```

Finally, we can also output our function to an HDF5 file, similar to the built-in functions to output selected field components, and so on. The following outputs an HDF5 file consisting of our function `f` evaluated at every point in the computational cell:

*Python*
```py
meep.Simulation.output_field_function("weird-function", [meep.Ex, meep.Hz, meep.Dielectric], f)
```

*Scheme*
```scm
(output-field-function "weird-function" (list Ex Hz Dielectric) f)
```

Here, the first argument is used for the name of the dataset within the HDF5, and is also used for the name of the HDF5 file itself plus a `.h5` suffix and a time stamp, unless you have specified the output file via `to_appended` or `to-appended` or other means.

The above example calls the integration, maximum, and output routines only once, at the current time. Often, you will want to pass them to `meep.Simulation.run(..., until=...)` (Python) or `run-until` (Scheme) instead, using `at_every` or `at-every` to print or output at periodic time intervals. In Scheme, a common mistake is to do something like the following:

```scm
(run-until 200 (at-every 1 (output-field-function "weird-function" (list Ex Hz Dielectric) f)))
```

This is **wrong**, and will cause Meep to exit with a strange error message. The reason is that the step functions you pass to `run-until` must be *functions*. For example, if you call `(run-until 200 output-hfield)`, `output-hfield` is the name of a *function* which `run-until` will call to output the field. The incorrect code above, however, first *calls* the function `output-field-function` to output an HDF5 file, and then passes the *result* of this function to `run-until`. Instead, you must write a new function which you can pass to `run-until`, like the following:

```scm
(define (my-weird-output) (output-field-function "weird-function" (list Ex Hz Dielectric) f))
(run-until 200 (at-every 1 my-weird-output))
```

Here, we have defined a function `my-weird-output` of no arguments that, when called, outputs our function `f`. We then pass this function to `run-until`. In contrast, our incorrect code above corresponds to passing `(my-weird-output)`, the *result* of calling `my-weird-output`, to `run-until`.

As described in [Synchronizing the Magnetic and Electric Fields](Synchronizing_the_Magnetic_and_Electric_Fields.md), because this example function combines electric and magnetic fields, we may want to synchronize them in time in order to compute this function more accurately, by wrapping it with `synchronized-magnetic`:

```scm
(run-until 200 (synchronized-magnetic (at-every 1 my-weird-output)))
```

See also the section "Writing Your Own Step Functions" in the [Python](Python_User_Interface.md#writing-your-own-step-functions) or [Scheme](Scheme_User_Interface.md#writing-your-own-step-functions) interface.
