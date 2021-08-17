# pfinder
App that generates fractals of polynomial roots given user input.  Built using a plotly Dash web framework running on gunicorn with computation performed using background Redis processes, hosted by Heroku.

After a short server activation time (this app runs on a free server which sleeps if not in use), the following interface will appear:

![screenshot](/assets/pfinder_screenshot.png)

submitting an equation (real values only at present) produces a high-resolution (up to 2400 by 1600 pixel) image, for example the default equation 'x^7.14-x-1' results in

![cover](/assets/pfinder_cover.png)

Complex-valued constants and exponents are now supported as well, for example:
'x^(7.11+0.1i)-x-(1-0.2i)'
Complex numbers may be entered with or without real values (ie (2i) is accetable) but must exist in single parentheses.

To learn more, see [this page](https://blbadger.github.io/polynomial-roots.html).  If this sort of thing is of interest, see also the Jenerator (https://github.com/blbadger/jenerator)


