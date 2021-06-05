"""Input and output file/directory utility functions."""

import os


def get_data(filename: str) -> str:
    """Gets the relative file path for the provided file."""
    return "./data/" + filename


def make_output_dir():
    """Creates the output file directory, if it doesn't exist."""
    if not os.path.exists("./output"):
        os.mkdir("./output")


def get_output(filename: str) -> str:
    """Gets the relative file path for the provided output file."""
    return "./output/" + filename
