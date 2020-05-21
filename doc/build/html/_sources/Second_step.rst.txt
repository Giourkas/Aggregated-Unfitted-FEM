For the second step, after we have distincted which elements are inside, outside and cut, the algorithm includes information about the adjacency of every cut cell to its closest interior cell (ref. :py:mod:`cluster_info`) defined by a distance (=the smallest number of cells needed for every cut cell to reach one interior cell) and it creates a cluster by grouping all cells that were mapped to the same inner cell.
The following image shows exactly one cluster, where the blue one inside is the inner cell and the greens are the cells cut by the boundary and clustered with the blue:

|one_cluster|

.. |one_cluster| image:: clustering2.png
   :height: 200
   :width: 200

In detail, we are starting with the interior cells and iteratively adding adjacent cut cells to it. This continues until all cut cells are clustered with an interior (and possibly other cut cells). This process can be shown in the below four steps.
