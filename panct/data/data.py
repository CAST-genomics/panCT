from __future__ import annotations
import os
import gzip
from pathlib import Path
from collections import namedtuple
from abc import ABC, abstractmethod
from logging import getLogger, Logger
from typing import Iterator, IO, Any, Type


class Data(ABC):
    """
    Abstract class for accessing read-only data files

    Attributes
    ----------
    data : np.array
        The contents of the data file, once loaded
    log: Logger
        A logging instance for recording debug statements.
    """

    def __init__(self, log: Logger = None):
        self.data = None
        self.log = log or getLogger(self.__class__.__name__)
        super().__init__()

    def __repr__(self):
        return self.data.__repr__()

    @classmethod
    @abstractmethod
    def read(cls: Type[Data], fname: Path | str) -> Data:
        """
        Read the file contents and perform any recommended pre-processing

        Parameters
        ----------
        fname : Path | str
            The name of a file to load data from

        Returns
        -------
        An initialized instance of Data with the data from fname loaded
        """
        pass

    # @abstractmethod
    # def __iter__(self) -> Iterator[namedtuple]:
    #     """
    #     Return an iterator over the raw file contents

    #     Yields
    #     ------
    #     Iterator[namedtuple]
    #         An iterator over each line in the file, where each line is encoded as a
    #         namedtuple containing each of the class properties
    #     """
    #     pass

    @staticmethod
    def hook_compressed(filename: Path | str, mode: str) -> gzip.GzipFile | IO[Any]:
        """
        A utility to help open files regardless of their compression

        Based off of python's fileinput.hook_compressed and copied from
        https://stackoverflow.com/a/64106815/16815703

        Parameters
        ----------
        filename : Path | str
            The path to the file
        mode : str
            Either 'r' for read or 'w' for write

        Returns
        -------
        gzip.GzipFile | IO[Any]
            The resolved file object
        """
        if "b" not in mode:
            mode += "t"
        ext = os.path.splitext(filename)[1]
        if ext == ".gz":
            return gzip.open(filename, mode)
        else:
            return open(filename, mode)
