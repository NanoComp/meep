--8<-- "AdjointSolver/AdjointDocumentationStyleHeader.md"

---
# Design optimization with the <span class=SC>meep</span> adjoint solver: A reference manual
---

## Defining your problem: Creating a subclass of `OptimizationProblem`

In the [Overview](Overview.md) discussion we ticked off
the various data items needed to define a `meep.adjoint` problem.
You will communicate this information to the solver
by writing a python script implementing a subclass
of `OptimizationProblem.` This is a base class implemented
by `meep.adjoint` that knows how to do various general things
involving <span class=SC>meep</span> geometries and objective
functions, but which is lacking crucial information from you, 
without which it can't do anything on its own.
(That is to say, `OptimizationProblem` is an
[abstract base class](https://docs.python.org/3/library/abc.html)
coining two pure virtual methods that your derived class
must override.

### Mandatory and optional class-method overrides

More specifically, your subclass of `OptimizationProblem`
must furnish implementations of the following two pure virtual
methods left

Having

## Evaluating your objective function

## Optimizing your geometry
