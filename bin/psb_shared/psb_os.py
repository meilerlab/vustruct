#!/usr/bin/env python3
import contextlib
import logging
import os

"""
Define a set of extraordinarily useful OS wrappers
"""

# def init_psb_os(dict)

def makedirs(name,exist_ok = True):
    """
    os.makedirs(
    """
    pass

LOGGER = logging.getLogger(__name__)

@contextlib.contextmanager
def pushdir(new_dir: str):
    """
    @new_dir: Directory to change to for duraiton of hte context
    Create a context where new_dir is the current directory until the context is edited.
    """
    previous_dir = os.getcwd()

    os.chdir(new_dir)

    LOGGER.info("chdir() from %s to %s", previous_dir,new_dir)

    try:
        yield
    finally:
        os.chdir(previous_dir)
        LOGGER.info("chdir() back to %s", previous_dir)
