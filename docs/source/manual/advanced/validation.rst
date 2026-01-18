Validation
===========

PhyloZoo includes a comprehensive validation system that ensures objects always represent 
valid phylogenetic structures. Validation occurs automatically upon object creation and 
can be controlled through decorators and context managers.

Automatic Validation
--------------------

By default, network objects are validated upon creation:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # This will raise PhyloZooNetworkStructureError due to cycle
   try:
       network = DirectedPhyNetwork(edges=[(1, 2), (2, 1)])
   except PhyloZooNetworkStructureError:
       print("Invalid network structure detected")

Validation ensures that:
* Network structure is valid (no cycles, proper connectivity)
* Node degrees satisfy constraints
* Edge attributes are valid (branch lengths, gamma values, etc.)
* Network properties are consistent

Validation Decorators
---------------------

The validation system uses decorators to control when validation occurs:

.. code-block:: python

   from phylozoo.utils.validation import validation_aware
   
   @validation_aware(allowed=["validate"], default=["validate"])
   class MyNetwork:
       def validate(self):
           # Validation logic
           pass

The ``@validation_aware`` decorator specifies which methods trigger validation and which 
methods are allowed to skip validation.

Disabling Validation
--------------------

For performance-critical code, validation can be temporarily disabled:

.. code-block:: python

   from phylozoo.utils.validation import no_validation
   
   with no_validation():
       # Operations that skip validation
       network = DirectedPhyNetwork(edges=[...])
       # ... perform operations ...

**Warning**: Only disable validation when you are certain the operations produce valid 
structures. Invalid networks can cause errors in downstream operations.

Custom Validation
-----------------

You can add custom validation logic to network classes:

.. code-block:: python

   class MyNetwork(DirectedPhyNetwork):
       def _validate_custom_constraint(self):
           # Custom validation logic
           if some_condition:
               raise PhyloZooNetworkStructureError("Custom constraint violated")

Custom validation methods should be named ``_validate_*`` and will be called automatically 
during validation.

Validation Methods
------------------

Network classes provide several validation methods:

**validate()**
   Performs full validation of the network structure, node degrees, and attributes. 
   Raises appropriate exceptions if validation fails.

**is_valid()**
   Checks if the network is valid without raising exceptions. Returns ``True`` if valid, 
   ``False`` otherwise.

**validate_structure()**
   Validates only the network structure (connectivity, cycles, etc.).

**validate_degrees()**
   Validates only node degree constraints.

**validate_attributes()**
   Validates only edge and node attributes.

Example: Custom Validation
---------------------------

.. code-block:: python

   from phylozoo import DirectedPhyNetwork, PhyloZooNetworkStructureError
   from phylozoo.utils.validation import validation_aware
   
   @validation_aware(allowed=["validate"], default=["validate"])
   class CustomNetwork(DirectedPhyNetwork):
       def _validate_max_level(self):
           """Ensure network level does not exceed maximum."""
           if self.level() > 3:
               raise PhyloZooNetworkStructureError(
                   f"Network level {self.level()} exceeds maximum of 3"
               )
       
       def validate(self):
           # Call parent validation
           super().validate()
           # Call custom validation
           self._validate_max_level()

Warnings
--------

PhyloZoo also provides custom warning classes for non-critical issues:

.. code-block:: python

   from phylozoo import (
       PhyloZooWarning,
       PhyloZooEmptyNetworkWarning,
       PhyloZooSingleNodeNetworkWarning,
   )

Warnings are issued for:
* Empty networks (technically valid but may not be useful)
* Single-node networks (valid but may not be useful for phylogenetic analysis)
* Other non-critical issues

Performance Considerations
---------------------------

**Validation Overhead**
   Validation adds computational overhead to object creation. For performance-critical 
   code, consider disabling validation when you're certain operations produce valid structures.

**Caching**
   PhyloZoo uses ``@cached_property`` for expensive computations. Once validated, network 
   properties are cached for subsequent access.

**Early Validation**
   Validation occurs early (upon object creation), catching errors before they propagate 
   to downstream operations. This makes debugging easier.

Best Practices
--------------

1. **Keep validation enabled**: Only disable validation when absolutely necessary for 
   performance. Validation catches errors early and makes debugging easier.

2. **Use custom validation**: Add custom validation methods for domain-specific constraints.

3. **Handle validation errors**: Always handle validation exceptions appropriately. Don't 
   suppress validation errors unless you're certain they're acceptable.

4. **Validate after modifications**: If you modify a network (by creating a new instance), 
   validation will occur automatically. For custom modifications, consider calling ``validate()`` 
   explicitly.

5. **Use is_valid() for checks**: When you need to check validity without raising exceptions, 
   use ``is_valid()`` instead of catching exceptions.

6. **Respect validation decorators**: When extending network classes, respect the validation 
   decorator settings and don't bypass validation unnecessarily.

Example: Validation Workflow
------------------------------

.. code-block:: python

   from phylozoo import DirectedPhyNetwork, PhyloZooNetworkError
   from phylozoo.utils.validation import no_validation
   
   # Normal creation with validation
   try:
       network = DirectedPhyNetwork(edges=[(1, 2), (2, 3)])
       print("Network is valid")
   except PhyloZooNetworkError as e:
       print(f"Validation failed: {e}")
   
   # Performance-critical code with validation disabled
   with no_validation():
       # Create many networks quickly (assumes they're valid)
       networks = [
           DirectedPhyNetwork(edges=[(i, i+1)]) 
           for i in range(100)
       ]
   
   # Validate after batch creation
   for network in networks:
       if not network.is_valid():
           print(f"Invalid network detected: {network}")

.. seealso::
   For more information on:
   * Exceptions: :doc:`Exceptions <../exceptions>`
   * Network structure: :doc:`Networks (Advanced) <../core/networks/advanced>`
   * Validation API: See ``phylozoo.utils.validation`` module documentation
