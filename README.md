## Pfinder: polynomial root finder

### How to use 

The Pfinder generates fractals with polynomial root-finding methods that color the complex plane depending on how long it takes for an initial point to find any root.  After following the link in the upper right hand corner of this page, a short server activation time ensues (this app currently runs on a free server which sleeps if not in use) and the following interface will appear:

![screenshot](/assets/pfinder_screenshot.png)

`Real Bounds` and `Imaginary Bounds` specify the location of the horizontal and vertical axes, respectively, in the complex plane: for example, imaginary bounds `-1, 1` results in an image generated for values between `-i` (at the bottom of the map) and `i` (top).  Inputs may or may not contain spaces, but must have a single comma `,` between boundary elements, and the first element must be less than the second.

The `Choose color map:` dropdown specify the colors assigned to the array of the number of iterations before arriving at a root.  The`Choose resolution:` dropdown menu specifies two integers, which are the number of values calculated between the `Real Bounds` and the `Imaginary Bounds`, respectively. There are approximately equal to the resolution of the image produced.

Any of Newton's method, Halley's method, or the Secant method (with the initial guess plotted and the second guess being half way between the origin `0 + 0i` and the first guess) can be selected in the `Choose method:` dropdown.

The `Maximum iterations` input selection field records the maximum number of iterations of the root method before the programs halts, and can be fed any value between `0` and `300`.  Roughly speaking, more iterations are required for higher resolution images.  Note that the different root methods require different iterations to make interesting maps, with the Secant method requiring the most and Halley's the fewest.

The last and most important input is the `Specify Equation` field, where the user specifies the equation that is then fed to the selected root method in order to generate a map of quickly- versus slowly-converging initial values. Inputs should contain no spaces, and only the characters `0123456789.ie^+-` should be entered.  The equation is parsed and rendered in MathJax to the right of the input field, and checking to make sure that the rendered equation matches the intended input is helpful to make sure that the parser recognizes the input properly.  

The example input equation `x^7.14-x-1` results in

![cover](/assets/pfinder_cover.png)

Complex-valued constants and exponents are now supported as well, for example:

`x^(7.11+0.1i)-x-(1-0.2i)`

Complex numbers may be entered with or without real values, ie `(2i)` is acceptable alternative to `(0+2i)`, but all complex or negative contants and exponents must exist in single parentheses.


### Behind the scenes

The Pfinder is a flask app using a Dash interface for front end layout and callbacks, and is styled in CSS.  Array-based computation necessary for image generation is handed off to a Redis server via a Redis Queue message broker, and the Redis server is pinged every two seconds by the app to see if computation is complete.  The actual computation is performed using Numpy, and the resulting array is transformed into an image via matplotlib and saved to a temporary memory buffer as a bytestring using a base64 binary encoding.  When the Redis job is complete, it is fetched and the bytestring is decoded into a PNG that can be opened in a separate page for maximum resolution.  

This system of background Redis processes circumvents the problem of long computation times faced at higher resolutions.  Heroku, Azure, and most other cloud PaaS providers have hard time limits (30s in this case) for safety and efficiency conerns, meaning that the long computations required to generate high-resolution images of complex-number arrays would simply time out without this system in place.  Timeouts do not occur even when computations run for more than 5 minutes with the current configuration because the app is continually 'active' as it pings the Redis server.  

The generation of a polynomial root map only occurs when the `CLICK TO RUN` button is pressed, but real-time callbacks occurs for rendering the specified equation with MathJax (with a delay of ~100ms).  This process does not employ redis but is computed in the flask app using helper functions, which allows for new equations to be rendered whilst a prior equation map is generated.


### Learn More

To learn more, see [this page](https://blbadger.github.io/polynomial-roots.html).  If this sort of thing is of interest to you, see also the [Jenerator](https://github.com/blbadger/jenerator)

