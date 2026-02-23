"""
Parallel execution utilities for computationally intensive operations.

This module provides a unified interface for parallelizing operations across
the PhyloZoo package. It supports multiple backends (sequential, threading,
multiprocessing) and can be used via function parameters.

Examples
--------
Basic usage with a function that accepts a ``parallel`` parameter:

    >>> from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
    >>>
    >>> # Use multiprocessing with 4 cores
    >>> result = some_parallel_function(
    ...     data,
    ...     parallel=ParallelConfig(
    ...         backend=ParallelBackend.MULTIPROCESSING,
    ...         n_jobs=4
    ...     )
    ... )
    >>>
    >>> # Use all available cores (auto-detect)
    >>> result = some_parallel_function(
    ...     data,
    ...     parallel=ParallelConfig(
    ...         backend=ParallelBackend.MULTIPROCESSING,
    ...         n_jobs=None  # or -1 for all cores
    ...     )
    ... )
    >>>
    >>> # Sequential execution (no parallelization)
    >>> result = some_parallel_function(
    ...     data,
    ...     parallel=ParallelConfig(backend=ParallelBackend.SEQUENTIAL)
    ... )

Using with combinations/iterables:

    >>> from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
    >>> import itertools
    >>> 
    >>> def process_quartet(indices):
    ...     i, j, k, l = indices
    ...     # Process quartet...
    ...     return result
    >>> 
    >>> combinations = list(itertools.combinations(range(20), 4))
    >>> config = ParallelConfig(
    ...     backend=ParallelBackend.MULTIPROCESSING,
    ...     n_jobs=4
    ... )
    >>> executor = config.get_executor()
    >>> results = list(executor.map(process_quartet, combinations))
"""

from __future__ import annotations

import itertools
import os
from enum import Enum
from typing import Any, Callable, Iterator, Protocol, TypeVar

from .exceptions import PhyloZooValueError

T = TypeVar('T')
R = TypeVar('R')


class ParallelBackend(Enum):
    """
    Available parallelization backends.
    
    Attributes
    ----------
    SEQUENTIAL : str
        No parallelization - executes sequentially (default).
        Use for debugging or when overhead outweighs benefits.
    THREADING : str
        Thread-based parallelization. Good for I/O-bound operations
        or when sharing memory is important. Limited by Python's GIL
        for CPU-bound tasks.
    MULTIPROCESSING : str
        Process-based parallelization. Best for CPU-bound tasks that
        don't require shared memory. Bypasses Python's GIL.
    """
    
    SEQUENTIAL = "sequential"
    THREADING = "threading"
    MULTIPROCESSING = "multiprocessing"


class ParallelExecutor(Protocol):
    """
    Protocol for parallel execution backends.
    
    All executors must implement `map` and `starmap` methods that
    apply a function to items in an iterable, potentially in parallel.
    """
    
    def map(
        self,
        func: Callable[[T], R],
        iterable: Iterator[T] | list[T],
        chunksize: int | None = None,
    ) -> Iterator[R]:
        """
        Apply function to each item in iterable, potentially in parallel.
        
        Parameters
        ----------
        func : Callable[[T], R]
            Function to apply to each item.
        iterable : Iterator[T] | list[T]
            Items to process.
        chunksize : int | None, optional
            Number of items to process per worker (backend-dependent).
            By default None (backend chooses).
        
        Yields
        ------
        R
            Results from applying func to each item.
        """
        ...
    
    def starmap(
        self,
        func: Callable[..., R],
        iterable: Iterator[tuple[Any, ...]] | list[tuple[Any, ...]],
        chunksize: int | None = None,
    ) -> Iterator[R]:
        """
        Apply function to unpacked arguments from iterable, potentially in parallel.
        
        Parameters
        ----------
        func : Callable[..., R]
            Function to apply. Will receive unpacked arguments from tuples.
        iterable : Iterator[tuple[Any, ...]] | list[tuple[Any, ...]]
            Tuples of arguments to unpack and pass to func.
        chunksize : int | None, optional
            Number of items to process per worker (backend-dependent).
            By default None (backend chooses).
        
        Yields
        ------
        R
            Results from applying func to unpacked arguments.
        """
        ...


class SequentialExecutor:
    """
    Sequential (no parallelization) executor.
    
    Executes operations in order, one at a time. Useful for debugging
    or when parallelization overhead is not worth it.
    """
    
    def map(
        self,
        func: Callable[[T], R],
        iterable: Iterator[T] | list[T],
        chunksize: int | None = None,
    ) -> Iterator[R]:
        """
        Apply function sequentially to each item.
        
        Parameters
        ----------
        func : Callable[[T], R]
            Function to apply.
        iterable : Iterator[T] | list[T]
            Items to process.
        chunksize : int | None, optional
            Ignored for sequential execution. By default None.
        
        Yields
        ------
        R
            Results from applying func to each item.
        """
        return map(func, iterable)
    
    def starmap(
        self,
        func: Callable[..., R],
        iterable: Iterator[tuple[Any, ...]] | list[tuple[Any, ...]],
        chunksize: int | None = None,
    ) -> Iterator[R]:
        """
        Apply function sequentially to unpacked arguments.
        
        Parameters
        ----------
        func : Callable[..., R]
            Function to apply.
        iterable : Iterator[tuple[Any, ...]] | list[tuple[Any, ...]]
            Tuples of arguments to unpack.
        chunksize : int | None, optional
            Ignored for sequential execution. By default None.
        
        Yields
        ------
        R
            Results from applying func to unpacked arguments.
        """
        return itertools.starmap(func, iterable)


class ThreadingExecutor:
    """
    Thread-based parallel executor.
    
    Uses Python's threading module. Good for I/O-bound operations or
    when sharing memory is important. Limited by Python's Global Interpreter
    Lock (GIL) for CPU-bound tasks.
    """
    
    def __init__(self, n_jobs: int | None = None) -> None:
        """
        Initialize threading executor.
        
        Parameters
        ----------
        n_jobs : int | None, optional
            Number of worker threads. If None or -1, uses 1 thread
            (threading is typically not beneficial for CPU-bound tasks).
            By default None.
        """
        from concurrent.futures import ThreadPoolExecutor
        
        if n_jobs is None or n_jobs == -1:
            n_jobs = 1  # Threading typically not beneficial for CPU-bound
        elif n_jobs <= 0:
            raise PhyloZooValueError(f"n_jobs must be positive, got {n_jobs}")
        
        self.n_jobs = n_jobs
        self._executor = ThreadPoolExecutor(max_workers=n_jobs)
    
    def map(
        self,
        func: Callable[[T], R],
        iterable: Iterator[T] | list[T],
        chunksize: int | None = None,
    ) -> Iterator[R]:
        """
        Apply function to items using thread pool.
        
        Parameters
        ----------
        func : Callable[[T], R]
            Function to apply.
        iterable : Iterator[T] | list[T]
            Items to process.
        chunksize : int | None, optional
            Ignored for threading executor. By default None.
        
        Yields
        ------
        R
            Results from applying func to each item.
        """
        return self._executor.map(func, iterable)
    
    def starmap(
        self,
        func: Callable[..., R],
        iterable: Iterator[tuple[Any, ...]] | list[tuple[Any, ...]],
        chunksize: int | None = None,
    ) -> Iterator[R]:
        """
        Apply function to unpacked arguments using thread pool.
        
        Parameters
        ----------
        func : Callable[..., R]
            Function to apply.
        iterable : Iterator[tuple[Any, ...]] | list[tuple[Any, ...]]
            Tuples of arguments to unpack.
        chunksize : int | None, optional
            Ignored for threading executor. By default None.
        
        Yields
        ------
        R
            Results from applying func to unpacked arguments.
        """
        return self._executor.map(lambda args: func(*args), iterable)
    
    def __del__(self) -> None:
        """Clean up thread pool executor."""
        if hasattr(self, '_executor'):
            self._executor.shutdown(wait=False)


class MultiprocessingExecutor:
    """
    Process-based parallel executor.
    
    Uses Python's multiprocessing module. Best for CPU-bound tasks as it
    bypasses Python's Global Interpreter Lock (GIL). Each worker is a
    separate process with its own memory space.
    """
    
    def __init__(self, n_jobs: int | None = None) -> None:
        """
        Initialize multiprocessing executor.
        
        Parameters
        ----------
        n_jobs : int | None, optional
            Number of worker processes. If None or -1, uses all available
            CPU cores (os.cpu_count()). Must be positive. By default None.
        
        Raises
        ------
        PhyloZooValueError
            If n_jobs is 0 or negative (except -1).
        """
        from multiprocessing import Pool
        
        if n_jobs is None or n_jobs == -1:
            n_jobs = os.cpu_count() or 1
        elif n_jobs <= 0:
            raise PhyloZooValueError(f"n_jobs must be positive, got {n_jobs}")
        
        self.n_jobs = n_jobs
        self._pool = Pool(processes=n_jobs)
    
    def map(
        self,
        func: Callable[[T], R],
        iterable: Iterator[T] | list[T],
        chunksize: int | None = None,
    ) -> Iterator[R]:
        """
        Apply function to items using process pool.
        
        Parameters
        ----------
        func : Callable[[T], R]
            Function to apply. Must be picklable for multiprocessing.
        iterable : Iterator[T] | list[T]
            Items to process. Items must be picklable.
        chunksize : int | None, optional
            Number of items to send to each worker at once. If None,
            multiprocessing chooses an appropriate value. By default None.
        
        Yields
        ------
        R
            Results from applying func to each item.
        """
        return iter(self._pool.map(func, iterable, chunksize=chunksize))
    
    def starmap(
        self,
        func: Callable[..., R],
        iterable: Iterator[tuple[Any, ...]] | list[tuple[Any, ...]],
        chunksize: int | None = None,
    ) -> Iterator[R]:
        """
        Apply function to unpacked arguments using process pool.
        
        Parameters
        ----------
        func : Callable[..., R]
            Function to apply. Must be picklable.
        iterable : Iterator[tuple[Any, ...]] | list[tuple[Any, ...]]
            Tuples of arguments to unpack. Must be picklable.
        chunksize : int | None, optional
            Number of items to send to each worker at once. By default None.
        
        Yields
        ------
        R
            Results from applying func to unpacked arguments.
        """
        return iter(self._pool.starmap(func, iterable, chunksize=chunksize))
    
    def __del__(self) -> None:
        """Clean up process pool."""
        if hasattr(self, '_pool'):
            self._pool.close()
            self._pool.join()


class ParallelConfig:
    """
    Configuration for parallel execution.
    
    This class encapsulates all settings needed for parallel execution,
    including backend selection and worker count. Use this as a function
    parameter (Pattern A) to enable parallelization in PhyloZoo functions.
    
    Parameters
    ----------
    backend : ParallelBackend | str, optional
        Parallelization backend to use. By default ParallelBackend.SEQUENTIAL.
    n_jobs : int | None, optional
        Number of workers. Interpretation depends on backend:
        - SEQUENTIAL: Ignored
        - THREADING: Number of threads (None/-1 = 1)
        - MULTIPROCESSING: Number of processes (None/-1 = all CPU cores)
        By default None.
    chunksize : int | None, optional
        Number of items to process per worker batch. Backend-dependent.
        By default None (backend chooses).
    
    Examples
    --------
    >>> from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
    >>> 
    >>> # Use 4 CPU cores
    >>> config = ParallelConfig(
    ...     backend=ParallelBackend.MULTIPROCESSING,
    ...     n_jobs=4
    ... )
    >>> 
    >>> # Use all available cores
    >>> config = ParallelConfig(
    ...     backend=ParallelBackend.MULTIPROCESSING,
    ...     n_jobs=None  # or -1
    ... )
    >>> 
    >>> # Sequential execution
    >>> config = ParallelConfig(backend=ParallelBackend.SEQUENTIAL)
    """
    
    def __init__(
        self,
        backend: ParallelBackend | str = ParallelBackend.SEQUENTIAL,
        n_jobs: int | None = None,
        chunksize: int | None = None,
    ) -> None:
        if isinstance(backend, str):
            try:
                self.backend = ParallelBackend(backend)
            except ValueError:
                raise PhyloZooValueError(
                    f"Unknown backend: {backend}. "
                    f"Must be one of: {', '.join(b.value for b in ParallelBackend)}"
                )
        else:
            self.backend = backend
        
        self.n_jobs = n_jobs
        self.chunksize = chunksize
    
    def get_executor(self) -> ParallelExecutor:
        """
        Get executor instance based on backend configuration.
        
        Returns
        -------
        ParallelExecutor
            Executor instance ready to use.
        
        Raises
        ------
        PhyloZooValueError
            If backend is not recognized.
        """
        if self.backend == ParallelBackend.SEQUENTIAL:
            return SequentialExecutor()
        elif self.backend == ParallelBackend.THREADING:
            return ThreadingExecutor(n_jobs=self.n_jobs)
        elif self.backend == ParallelBackend.MULTIPROCESSING:
            return MultiprocessingExecutor(n_jobs=self.n_jobs)
        else:
            raise PhyloZooValueError(f"Unknown backend: {self.backend}")
    
    def __repr__(self) -> str:
        """String representation of configuration."""
        return (
            f"ParallelConfig(backend={self.backend.value}, "
            f"n_jobs={self.n_jobs}, chunksize={self.chunksize})"
        )
