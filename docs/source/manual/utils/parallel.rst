Parallel Execution
==================

The :mod:`phylozoo.utils.parallel` module provides the standard interface for
parallel execution across PhyloZoo. It defines a small set of backends
(sequential, threading, multiprocessing) behind a single
:class:`~phylozoo.utils.parallel.ParallelConfig` configuration object and a
common executor protocol, so that any function that wants to support
parallelization exposes the same API to its users.

All classes on this page can be imported from the parallel module:

.. code-block:: python

   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend

.. note::
   No PhyloZoo function currently exposes a ``parallel`` parameter. This module
   is the intended way to introduce parallelization in future implementations or
   for dependent packages to utilize, so that all parallelized functions share a uniform interface.

How parallelization is organized
--------------------------------

Parallel execution in PhyloZoo follows a single pattern: a function that
supports parallel execution accepts a :class:`~phylozoo.utils.parallel.ParallelConfig`
(or ``None`` for the default sequential behavior), obtains an executor from it,
and uses the executor's :meth:`map` / :meth:`starmap` to dispatch work. Choosing
a backend, the number of workers, and the chunk size is therefore done in one
place (the config), and the function itself does not need to know about
threads, processes, or pools.

Backends
^^^^^^^^

Three backends are available through
:class:`~phylozoo.utils.parallel.ParallelBackend`:

- ``SEQUENTIAL`` — no parallelization; executes one item at a time. The default,
  and the right choice for debugging or when parallel overhead would outweigh
  the benefits.
- ``THREADING`` — uses a thread pool (:class:`concurrent.futures.ThreadPoolExecutor`).
  Good for I/O-bound work or when worker memory must be shared; limited by
  Python's GIL for CPU-bound tasks.
- ``MULTIPROCESSING`` — uses a process pool (:class:`multiprocessing.Pool`).
  Best for CPU-bound tasks; bypasses the GIL but requires that the function
  and its arguments be picklable.

Each backend is implemented as a small executor class
(:class:`~phylozoo.utils.parallel.SequentialExecutor`,
:class:`~phylozoo.utils.parallel.ThreadingExecutor`,
:class:`~phylozoo.utils.parallel.MultiprocessingExecutor`) that conforms to the
:class:`~phylozoo.utils.parallel.ParallelExecutor` protocol with ``map`` and
``starmap`` methods.

Configuring parallel execution
------------------------------

Construct a :class:`~phylozoo.utils.parallel.ParallelConfig` to describe how
work should be dispatched.

**Sequential (default)**

.. code-block:: python

   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend

   config = ParallelConfig(backend=ParallelBackend.SEQUENTIAL)

**Multiprocessing with a fixed number of workers**

.. code-block:: python

   config = ParallelConfig(
       backend=ParallelBackend.MULTIPROCESSING,
       n_jobs=4,
   )

**Multiprocessing using all available cores**

Pass ``n_jobs=None`` (or ``-1``) to use every available CPU core:

.. code-block:: python

   config = ParallelConfig(
       backend=ParallelBackend.MULTIPROCESSING,
       n_jobs=None,
   )

**Threading**

.. code-block:: python

   config = ParallelConfig(
       backend=ParallelBackend.THREADING,
       n_jobs=4,
   )

The backend can also be supplied as a string (``"sequential"``,
``"threading"``, ``"multiprocessing"``) for convenience.

Using an executor directly
--------------------------

Once a config is built, call
:meth:`~phylozoo.utils.parallel.ParallelConfig.get_executor` to obtain an
executor and use its :meth:`map` or :meth:`starmap` methods:

.. code-block:: python

   import itertools
   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend

   def process_quartet(indices: tuple[int, int, int, int]) -> float:
       i, j, k, l = indices
       # ... do work ...
       return 0.0

   combinations = list(itertools.combinations(range(20), 4))

   config = ParallelConfig(
       backend=ParallelBackend.MULTIPROCESSING,
       n_jobs=4,
   )
   executor = config.get_executor()
   results = list(executor.map(process_quartet, combinations))

For functions that take multiple positional arguments, use
:meth:`~phylozoo.utils.parallel.ParallelExecutor.starmap` with an iterable of
argument tuples:

.. code-block:: python

   def pair_distance(a: int, b: int) -> float:
       # ... do work ...
       return 0.0

   pairs = [(0, 1), (0, 2), (1, 2)]
   results = list(executor.starmap(pair_distance, pairs))

Pattern for PhyloZoo functions
------------------------------

When a PhyloZoo function wants to support parallel execution, it should accept
a ``parallel`` keyword argument typed as
:class:`~phylozoo.utils.parallel.ParallelConfig` ``| None``, default to a
sequential config when ``None`` is given, and dispatch its inner loop through
the executor:

.. code-block:: python

   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend

   def my_parallel_function(
       data: list[int],
       parallel: ParallelConfig | None = None,
   ) -> list[float]:
       if parallel is None:
           parallel = ParallelConfig(backend=ParallelBackend.SEQUENTIAL)
       executor = parallel.get_executor()
       return list(executor.map(_process_item, data))

This keeps the choice of backend and worker count in the caller's hands while
the function body stays agnostic to the parallelization strategy.

.. warning::
   For the multiprocessing backend, both ``func`` and every item in
   ``iterable`` must be picklable. Module-level functions are picklable;
   closures and lambdas are not.

See Also
--------

- :doc:`Parallel API <../../api/utils/parallel>` — Full class and method reference.
- :doc:`Exceptions <exceptions>` — :class:`~phylozoo.utils.exceptions.general.PhyloZooValueError` is raised for invalid backend names or non-positive worker counts.
