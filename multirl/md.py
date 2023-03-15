"""Functions related to running and analyzing molecular dynamic trajectories"""
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path

from multirl.models import Sequence
from multirl.utils import pin_then_run
from multirl.sim.run import sim_eval


class MolecularDynamicsResult(str):
    """Result from a molecular dynamics calculation (placeholder)"""


def batch_run_molecular_dynamics(sequences: list[Sequence], md_yml: Path) -> list[MolecularDynamicsResult]:
    """Run molecular dynamics on each member of a batch of sequences

    Args:
        sequences: Batch of sequences to be evaluated. Will run each of them in parallel with each other
    Returns:
        Scores for each
    """

    # Set the number of workers based on the batch size
    n_ranks = len(sequences)

    with ProcessPoolExecutor(n_ranks) as exc:
        # Submit all sequences
        futures = [
            exc.submit(pin_then_run, run_molecular_dynamics, rank, n_ranks, seq, md_yml)
            for rank, seq in enumerate(sequences)
        ]

        # Return all results
        return [f.result() for f in futures]


def run_molecular_dynamics(sequence: Sequence, md_yml: Path) -> MolecularDynamicsResult:
    """Perform a molecular dynamics calculation on a sequence

    Args:
        sequence: Sequence to be evaluated
        md_yml: Yaml file that stores MD simulation setup
    Returns:
        Result of the molecular dynamics
    """
    assert sequence is not None
    pdb = ESMFold(sequence)
    rmsf = sim_eval(pdb, md_yml)
    return MolecularDynamicsResult(os.environ.get("CUDA_VISIBLE_DEVICES"))
