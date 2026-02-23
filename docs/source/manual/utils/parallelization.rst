Parallelization
===============

PhyloZoo includes a unified parallelization interface for computationally intensive operations. 
Functions that accept a ``parallel`` parameter can utilize multiple CPU cores.

Basic Usage
-----------

The parallelization system uses a function parameter pattern. Functions that support 
parallelization accept an optional ``parallel`` parameter:

.. code-block:: python

   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
   
   # Use multiprocessing with 4 cores
   result = some_parallel_function(
       data,
       parallel=ParallelConfig(
           backend=ParallelBackend.MULTIPROCESSING,
           n_jobs=4
       )
   )
   
   # Use all available cores (auto-detect)
   result = some_parallel_function(
       data,
       parallel=ParallelConfig(
           backend=ParallelBackend.MULTIPROCESSING,
           n_jobs=None  # or -1 for all cores
       )
   )
   
   # Sequential execution (no parallelization, default)
   result = some_parallel_function(
       data,
       parallel=ParallelConfig(backend=ParallelBackend.SEQUENTIAL)
   )

Available Backends
------------------

**SEQUENTIAL** (default)
   No parallelization. Executes operations sequentially, one at a time. Useful for 
   debugging or when parallelization overhead outweighs benefits.

**THREADING**
   Thread-based parallelization. Good for I/O-bound operations or when sharing memory 
   is important. Limited by Python's Global Interpreter Lock (GIL) for CPU-bound tasks.

**MULTIPROCESSING**
   Process-based parallelization. Best for CPU-bound tasks as it bypasses Python's GIL. 
   Each worker is a separate process with its own memory space. Recommended for most 
   computationally intensive operations.

Worker Count (n_jobs)
----------------------

The ``n_jobs`` parameter controls the number of workers:

* ``None`` or ``-1``: Auto-detect (uses all available CPU cores for multiprocessing, 1 for threading)
* Positive integer: Use that many workers
* Must be positive (except -1 for "all cores")

.. code-block:: python

   import os
   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
   
   # Use all cores
   config = ParallelConfig(
       backend=ParallelBackend.MULTIPROCESSING,
       n_jobs=None  # or -1
   )
   
   # Use specific number of cores
   config = ParallelConfig(
       backend=ParallelBackend.MULTIPROCESSING,
       n_jobs=4
   )
   
   # Check available cores
   print(f"Available CPU cores: {os.cpu_count()}")

When to Use Parallelization
----------------------------

Parallelization is most beneficial for:

* **Embarrassingly parallel operations**: Operations that can be split into independent tasks
* **CPU-bound computations**: Computations limited by CPU speed, not I/O
* **Large datasets**: When processing many items (e.g., many quartets, many trees)

Parallelization may not help (or may even slow down) for:

* **Small datasets**: Overhead of process/thread creation exceeds computation time
* **I/O-bound operations**: Operations limited by disk/network speed
* **Operations with shared state**: When tasks need to share mutable data

Example: Processing items in parallel
-------------------------------------

.. code-block:: python

   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
   
   # Process multiple items in parallel
   items = [item1, item2, item3, item4]
   
   config = ParallelConfig(
       backend=ParallelBackend.MULTIPROCESSING,
       n_jobs=4
   )
   executor = config.get_executor()
   
   results = list(executor.map(process_item, items))

Best Practices
--------------

1. **Test performance**: Measure whether parallelization actually speeds up your workflow. 
   For small datasets, sequential execution may be faster due to overhead.

2. **Choose the right backend**: Use ``MULTIPROCESSING`` for CPU-bound tasks, ``THREADING`` 
   for I/O-bound tasks, and ``SEQUENTIAL`` for debugging or when overhead is too high.

3. **Use appropriate worker count**: Don't use more workers than CPU cores. Using more 
   workers than cores can actually slow down execution due to context switching overhead.

4. **Handle serialization**: Be aware that multiprocessing requires objects to be 
   serializable (picklable). PhyloZoo handles this automatically for supported operations.

5. **Consider memory usage**: Each worker process has its own memory space. For large 
   objects, this can significantly increase memory usage.

6. **Use sequential for debugging**: When debugging, use sequential execution to avoid 
   complications from parallel execution.

.. seealso::
   Parallelization API: See ``phylozoo.utils.parallel`` module documentation
