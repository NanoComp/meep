from typing import Dict, List, Optional

import numpy as np

import meep as mp

# Codes for different Meep time sinks, used by `mp.Simulation.time_spent_on()`.
# See
# https://meep.readthedocs.io/en/latest/Python_User_Interface/#simulation-time
# for more information.
TIMING_MEASUREMENT_IDS = {
    "connecting_chunks": mp.Connecting,
    "time_stepping": mp.Stepping,
    "boundaries_copying": mp.Boundaries,
    "mpi_all_to_all": mp.MpiAllTime,
    "mpi_one_to_one": mp.MpiOneTime,
    "field_output": mp.FieldOutput,
    "fourier_transform": mp.FourierTransforming,
    "mpb": mp.MPBTime,
    "near_to_farfield_transform": mp.GetFarfieldsTime,
    "other": mp.Other,
    "field_update_b": mp.FieldUpdateB,
    "field_update_h": mp.FieldUpdateH,
    "field_update_d": mp.FieldUpdateD,
    "field_update_e": mp.FieldUpdateE,
    "boundary_stepping_b": mp.BoundarySteppingB,
    "boundary_stepping_wh": mp.BoundarySteppingWH,
    "boundary_stepping_ph": mp.BoundarySteppingPH,
    "boundary_stepping_h": mp.BoundarySteppingH,
    "boundary_stepping_d": mp.BoundarySteppingD,
    "boundary_stepping_we": mp.BoundarySteppingWE,
    "boundary_stepping_pe": mp.BoundarySteppingPE,
    "boundary_stepping_e": mp.BoundarySteppingE,
}


def get_current_time_step(sim: mp.Simulation) -> int:
    """Gets the current time step from a Meep simulation object."""
    return sim.fields.t


class MeepTimingMeasurements:
    """Container for Meep timing measurements and performance information.

    See
    https://meep.readthedocs.io/en/latest/Python_User_Interface/#simulation-time
    for more information on the Meep timing measurements.

    Attributes:
      measurements: a dictionary of timing measurements, where each entry is a
        list containing timing measurements for each process.
      elapsed_time: the simulation's elapsed wall time.
      num_time_steps: the total number of time steps taken in the simulation.
      time_per_step: the wall time taken per step. Each element of this list is
        expected to correspond to a unique time step, but could also be averaged
        over several time steps.
      dft_relative_change: the relative change in the DFT fields, measured at each
        time step.
      overlap_relative_change: the relative change in the mode overlaps, measured
        at each time step.
      relative_energy: the relative energy in the simulation domain, measured
        at each time step.
      measurement_names: a list of measurement names stored in `measurements`.
      comm_efficiency: the communication efficiency, defined as the mean
        MPI/synchronization time divided by the sum of the mean time taken for
        time stepping and the mean time taken for fourier transforming. This is
        zero if no simulation work has been performed.
      comm_efficiency_one_to_one: MPI one-to-one communication efficiency, defined
        as the mean MPI one-to-one communication time divided by the sum of the
        mean time taken for time stepping and the mean time taken for fourier
        transforming. This is zero if no simulation work has been performed.
      comm_efficiency_all_to_all: MPI all-to-all communication efficiency, defined
        as the mean MPI all-to-all communication time divided by the sum of the
        mean time taken for time stepping and the mean time taken for fourier
        transforming. This is zero if no simulation work has been performed.
    """

    def __init__(
        self,
        measurements: Dict[str, List[float]],
        elapsed_time: float,
        num_time_steps: int,
        time_per_step: List[float],
        dft_relative_change: List[float],
        overlap_relative_change: List[float],
        relative_energy: List[float],
    ):
        self.measurements = measurements
        self.elapsed_time = elapsed_time
        self.num_time_steps = num_time_steps
        self.time_per_step = time_per_step
        self.dft_relative_change = dft_relative_change
        self.overlap_relative_change = overlap_relative_change
        self.relative_energy = relative_energy

    @classmethod
    def new_from_simulation(
        cls,
        sim: mp.Simulation,
        elapsed_time: Optional[float] = -1,
        time_per_step: Optional[List[float]] = None,
        dft_relative_change: Optional[List[float]] = None,
        overlap_relative_change: Optional[List[float]] = None,
        relative_energy: Optional[List[float]] = None,
    ) -> "MeepTimingMeasurements":
        """Creates a new `MeepTimingMeasurements` from a Meep simulation object.

        Usage example:
        ```
        sim = mp.Simulation(...)
        sim.init_sim()
        start_time = time.time()
        sim.run(...)
        measurements = meep_timing.MeepTimingMeasurements.new_from_simulation(
            sim=sim,
            elapsed_time=time.time() - start_time,
        )
        ```

        Args:
          sim: a Meep simulation object that has previously been run.
          elapsed_time: the wall time that has elapsed while running the simulation.
          time_per_step: the wall time taken per step. Each element of this list is
            expected to correspond to a unique time step, but could also be averaged
            over several time steps.
          dft_relative_change: the relative change in the DFT fields, measured at
            each time step.
          overlap_relative_change: the relative change in the mode overlaps,
            measured at each time step.
          relative_energy: the relative energy in the simulation domain, measured
            at each time step.

        Returns:
          the resulting `MeepTimingMeasurements`
        """
        measurements = {
            name: sim.time_spent_on(timing_id).tolist()
            for name, timing_id in TIMING_MEASUREMENT_IDS.items()
        }

        return cls(
            measurements=measurements,
            elapsed_time=elapsed_time,
            num_time_steps=get_current_time_step(sim),
            time_per_step=time_per_step if time_per_step is not None else [],
            dft_relative_change=dft_relative_change
            if dft_relative_change is not None
            else [],
            overlap_relative_change=overlap_relative_change
            if overlap_relative_change is not None
            else [],
            relative_energy=relative_energy if relative_energy is not None else [],
        )

    @property
    def measurement_names(self) -> List[str]:
        return list(self.measurements)

    @property
    def comm_efficiency(self) -> float:
        computation_time = np.mean(self.measurements["time_stepping"]) + np.mean(
            self.measurements["fourier_transform"]
        )
        if computation_time == 0:
            return 0
        else:
            return (
                np.mean(
                    self.measurements["mpi_all_to_all"]
                    + self.measurements["mpi_one_to_one"]
                )
                / computation_time
            )

    @property
    def comm_efficiency_one_to_one(self) -> float:
        computation_time = np.mean(self.measurements["time_stepping"]) + np.mean(
            self.measurements["fourier_transform"]
        )
        if computation_time == 0:
            return 0
        else:
            return np.mean(self.measurements["mpi_one_to_one"]) / computation_time

    @property
    def comm_efficiency_all_to_all(self) -> float:
        computation_time = np.mean(self.measurements["time_stepping"]) + np.mean(
            self.measurements["fourier_transform"]
        )
        if computation_time == 0:
            return 0
        else:
            return np.mean(self.measurements["mpi_all_to_all"]) / computation_time
