---
title: Synchronizing the magnetic and electric fields
permalink: /Synchronizing_the_magnetic_and_electric_fields/
---

In the finite-difference time-domain method, the electric and magnetic fields are stored at *different times* (and different positions in space), in a "[leap-frog](/w:Leapfrog_integration "wikilink")" fashion. At any given time-step \(t\) during the simulation, the **E** and **D** fields are stored at time \(t\), but the **H** and **B** fields are stored at time \(t-\Delta t/2\) (where \(\Delta t\) is the time-step size).

This means that when you output the electric and magnetic fields from a given time step, for example, the fields actually correspond to times \(\Delta t/2\) apart. For most purposes, this slight difference in time doesn't actually matter much, but it makes a difference when you compute quantities like the Poynting flux \(\mathbf{E}\times\mathbf{H}\) that combine electric and magnetic fields together, e.g. for the `output-poynting` function. If what you really want is the Poynting flux \(\mathbf{S}(t)\) at time *t*, then computing \(\mathbf{E}(t)\times\mathbf{H}(t-\Delta t/2)\) is slightly off from this — the error is of order \(O(\Delta t)\), or first-order accuracy. This is unfortunate, because the underlying FDTD method ideally can have second-order accuracy.

To improve the accuracy for computations involving both electric and magnetic fields, Meep provides a facility to synchronize the **H** and **B** fields with the **E** and **D** fields in time. Technically, what it does is to compute the magnetic fields at time \(t+\Delta t/2\) by performing part of a timestep, and then averaging those fields with the fields at time \(t-\Delta t/2\). This produces the magnetic fields at time *t* to second-order accuracy \(O(\Delta t^2)\), which is the best we can do in second-order FDTD. Meep also saves a copy of the magnetic fields at \(t-\Delta t/2\), so that it can restore those fields for subsequent timestepping.

Synchronization functions
-------------------------

All of this process is handled for you in Meep by a single step function: `synchronized-magnetic`. By wrapping this around your step functions, it ensures that those step functions are called with synchronized electric and magnetic fields (to second-order accuracy), while restoring the magnetic fields automatically for subsequent timestepping.

For example, if you do:

`(run-until 200 output-poynting output-tot-pwr)`

it outputs the Poynting vector and the total energy density in the electric and magnetic fields at each timestep, but it only does so to first-order accuracy because those computations combine unsynchronized electric and magnetic fields. Instead, if you do

`(run-until 200 (synchronized-magnetic output-poynting output-tot-pwr))`

it will output the same quantities, but more accurately because the fields will be synchronized. (Of course, **there is a price**: synchronizing the fields takes time, and also increases the memory usage in order to backup the unsynchronized fields.)

Alternatively, if you want to synchronize the magnetic and electric fields in some context other than that of a step function, e.g. you are doing some computation like `integrate-field-function` outside of the timestepping, you can instead call two lower-level functions. Before doing your computations, you should call `(meep-fields-synchronize-magnetic-fields` `fields)` to synchronize the magnetic fields with the electric fields, and after your computation you should call `(meep-fields-restore-magnetic-fields` `fields)` to restore the fields to their unsynchronized state for timestepping. (In the C++ interface, these correspond to `fields::synchronize_magnetic_fields` and `fields::restore_magnetic_fields`.) If you *don't* call `meep-fields-restore-magnetic-fields` before timestepping, then the fields will be re-synchronized after *every* timestep, which will greatly increase the cost of timestepping.

**Note**: in future versions of Meep, we may decide to synchronize the fields automatically whenever you output something like the Poynting vector or do another field computation that involves both magnetic and electric fields, but currently you must do this manually. (In any case, Meep does no additional work when you nest synchronization calls, so it is harmless to insert redundant field synchronizations.) The `flux-in-box` and `field-energy-in-box` routines are already automatically synchronized, however.

[Category:Meep](/Category:Meep "wikilink")