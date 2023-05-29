Examples
========

We shared two examples of the code.  

ABC grain 
---------

.. mat:function:: example1_ABCToyWithStepMobil()

This code runs grain growth simulation with simple anisotropic grain boundary characters. In particular, grains are grouped into 3 types, A,B,C. Grain boundaries between different groups have larger grain boundary energy and higher mobilities, while the same group has lower energy and mobiltiy.

.. mat:function:: example2_HAGB_poly()

This code runs grain growth simulation with [110] symmetric tilt grain boundary energies. Orientation values of grains are sampled from
sub-domain region with misorientation angle [0,70.6] degree from a reference state. Also, mobilities for small angle grain boundaries are set to lower than that of high angle grain boundaries. 

