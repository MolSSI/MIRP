.. _theory_gpt:

===========================
Gaussian Product Theorem
===========================

.. math::

   \gamma &= \alpha_1 + \alpha_2 \\
   \overline{AB}^2 &= (A_x - B_x)^2 + (A_y - B_y)^2 + (A_z - B_z)^2 \\
   P_x &= \frac{(\alpha_1 A_x + \alpha_2 B_x)}{\gamma} \\
   P_y &= \frac{(\alpha_1 A_y + \alpha_2 B_y)}{\gamma} \\
   P_z &= \frac{(\alpha_1 A_z + \alpha_2 B_z)}{\gamma} \\
   \overline{PA}_x &= P_x - A_x \qquad\qquad
   \overline{PA}_y = P_y - A_y \qquad\qquad
   \overline{PA}_z = P_z - A_z \\
   \overline{PB}_x &= P_x - B_x \qquad\qquad
   \overline{PB}_y = P_y - B_y \qquad\qquad
   \overline{PB}_z = P_z - B_z

   

---------
Functions
---------

.. doxygenfunction:: mirp_gpt
.. doxygenfunction:: mirp_gpt_mp
.. doxygenfunction:: mirp_gpt_interval

