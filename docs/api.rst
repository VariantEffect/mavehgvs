.. _api-docs:

mavehgvs API documentation
==========================

Variant objects
---------------

Each variant can be parsed into a variant object, which populates and exposes named
fields for each piece of the variant string.

.. autoclass:: mavehgvs.position.VariantPosition
   :members:

.. autoclass:: mavehgvs.variant.Variant
   :members:

Utility functions for regular expression patterns
-------------------------------------------------

.. automodule:: mavehgvs.patterns.util
   :members:

DNA pattern strings
-------------------

.. automodule:: mavehgvs.patterns.dna
   :members:

RNA pattern strings
-------------------

.. automodule:: mavehgvs.patterns.rna
   :members:

Protein pattern strings
-----------------------

.. automodule:: mavehgvs.patterns.protein
   :members:
