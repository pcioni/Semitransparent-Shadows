HOMEWORK 2: CLOTH & FLUID SIMULATION

NAME:  Anthony LaRosa



ESTIMATE OF # OF HOURS SPENT ON THIS ASSIGNMENT:  < 13 hrs >



COLLABORATORS AND OTHER RESOURCES: List the names of everyone you
talked to about this assignment and all of the resources (books,
online reference material, etc.) you consulted in completing this
assignment.

Phil Cioni
Jeramey Tyler

Remember: Your implementation for this assignment must be done on your
own, as described in "Academic Integrity for Homework" handout.



DISCUSS STABILITY OF YOUR CLOTH SIMULATION:
Very stable, works with all test cases, including my "custom.txt" cloth example.
Usually timesteps must be in thousandths (i.e. 0.001, 0.005)


DESCRIBE YOUR NEW CLOTH TEST SCENE &
THE COMMAND LINE TO RUN YOUR EXAMPLE:

It is an extension of table cloth, however the anchors are at different elevations.
It shows an interesting ridge across two diagonals, where it should rest in a way
where the cloth's springs are underextended.

./simulation -cloth ../src/custom.txt -timestep 0.001


DISCUSS THE ACCURACY & STABILITY OF YOUR FLUID SIMULATION:
Pretty accurate, but there seems to be an off by half error that I am unable to track down.
Despite this, it works in 3 dimensions.


KNOWN BUGS IN YOUR CODE:
Please be concise!
Timesteps in the hundreths break my cloth.
The "zero" boundaries of the fluid are not functional due to suspected off by half error.
Incompressiblity explodes, and is skipped.

NEW FEATURES OR EXTENSIONS FOR EXTRA CREDIT:
Include instructions for use and test cases and sample output as appropriate.
Fluid works in 3D.
Midpoint Integration and RungeKutta integration is implemented, as well as command line functionality.
Use -method argument, followed by one of the following arguments (defaults to euler):
  euler
  midpoint
  rungekutta

The terminal will display a line of output stating which method it is using. If it says Euler
when you typed something else, you may have mistyped the method.

./simulation -cloth ../src/small_cloth.txt -timestep 0.001 -method rungekutta
