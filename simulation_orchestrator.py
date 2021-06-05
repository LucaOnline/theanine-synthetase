from subprocess import Popen, DEVNULL

from options import SIMULATION_COUNT


def orchestrate_simulations(parallelism: int = 1):
    if SIMULATION_COUNT % parallelism != 0:
        print("Simulations cannot be divided evenly into requested instance count!")
        exit(1)

    simulations = [
        Popen(
            [
                "python",
                "simulation.py",
                "--id",
                str(i),
                "--trials",
                str(SIMULATION_COUNT // 8),
            ],
            stdout=DEVNULL,
        )
        for i in range(8)
    ]

    for simulation in simulations:
        simulation.wait()


if __name__ == "__main__":
    orchestrate_simulations(parallelism=8)
