The user interface code for the PSB Pipeline is taken from the NGL Webviewer
application.  The most changes from the NGL repository are reflected in the new
file examples/js/psb_gui.js which is used in place of the provided gui.js

The color_chains_consecutively module has routines to color chains not only consecutively
but also for alpha fold and pathprox colorings

Note that you cannot upgrade ngl.js without breaking chain coloring.  I have zero idea why

