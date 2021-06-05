# CSS383-Project1
Repository for CSS 383 Project 1

## Usage
Please refer to the [Makefile](./Makefile) for script usage.

`make install` installs the project dependencies from PyPi. You may want to create a [virtual environment](https://docs.python.org/3/library/venv.html) before running this command.

`make data` downloads and unpacks the datasets used in this project.

`make analyze` runs analyses on theanine synthetase and glutamine synthetase.

`make simulate` runs multiple simulation instances in parallel to determine how plausible it is that
the mutations between CsTSI and CsGSI occurred by random chance.

## Contributing
This project uses [Black](https://black.readthedocs.io/en/stable/) for Python code formatting.
