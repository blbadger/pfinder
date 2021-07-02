# app.py
# Responsive dashboard for exploring Julia sets

# import standard libraries
import io
import base64

# import third-party libraries
import numpy as np 
import matplotlib.pyplot as plt 
import dash
import dash_core_components as dcc 
import dash_html_components as html 
from dash.dependencies import Input, Output
import flask
import json
from rq.serializers import JSONSerializer
from dash.exceptions import PreventUpdate

from redis import Redis
from rq import Queue
from CalculateFaster import OptiCalculate
from rootmethod import newton, halley, secant

# import connection from worker.py
from worker import conn 
import time

q = Queue(connection=conn)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

MATHJAX_CDN = '''
https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/
MathJax.js?config=TeX-MML-AM_CHTML
'''

external_scripts = [
                    {'type': 'text/javascript',
                     'id': 'MathJax-script',
                     'src': MATHJAX_CDN,
                     },
                    ]


server = flask.Flask(__name__)
app = dash.Dash(server=server, 
				title='Polynomial roots', 
				external_stylesheets=external_stylesheets, 
				external_scripts=external_scripts)

job = ''

colors = {
	'background': '#fffff',
	'text': '#10110'
}

colormaps = [
			'viridis', 'plasma', 'inferno', 'magma', 'cividis',
			'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
			'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
			'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
			'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
			'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 
			'seismic','binary', 'gist_yarg', 'gist_gray', 'gray', 
			'bone', 'pink','spring', 'summer', 'autumn', 'winter', 
			'cool', 'Wistia','hot', 'afmhot', 'gist_heat', 'copper',
			'twilight', 'twilight_shifted', 'hsv','Pastel1', 'Pastel2', 
			'Paired', 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3',
			'tab10', 'tab20', 'tab20b', 'tab20c','flag', 'prism', 
			'ocean', 'gist_earth', 'terrain', 'gist_stern',
			'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
			'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral',
			'gist_ncar'
			]

methods = ["Newton's", "Halley's", "Secant"]

# page layout and inputs specified
app.layout = html.Div(
	style=
		{'backgroundColor': colors['background'], 
		'font-family': 'sans-serif'}, 
	children=[
	html.H1(
		children='Pfinder',
		style={
			'textAlign': 'center',
			'color': colors['text'],
			'margin-bottom': '1vh',
			'margin-top': '1vh'
		}
	),


	html.Div(
		children=[
		html.H4(
			children='Polynomial root finder',
			style={
				'textAlign': 'center',
				'color': colors['text'],
				'margin-bottom': '0vh',
				'margin-top': '0vh',
				'margin-left': '55vw',
				'width': '40vw',
				'height': '4vh',
				'display': 'inline-block',
				'vertical-align': 'top'
			}
		),
		html.Br(),

		html.Div(
			children=[
				html.Label('Real bounds'),
					dcc.Input(
					id='rbounds',
					type='text',
					value='-1.2, 1.2',
					style={'margin-top': '1vh',
							'width': '14vw'})
			], 
			style={
			'display':'inline-block',
			'width': '14vw',
			'margin-left': '2vw', 
			'margin-right': '2vw',
			'margin-top': '0vh',
			'text-align': 'center',
			'vertical-align': 'bottom'
		}),

		html.Div(
			children=[
				html.Label('Imaginary bounds'),
					dcc.Input(
					id='ibounds',
					type='text',
					value='-1, 1',
					style={'margin-top': '1vh',
							'width': '14vw'})
			], 
			style={
			'display':'inline-block',
			'width': '14vw',
			'margin-left': '1vw', 
			'margin-right': '1vw',
			'margin-top': '0vh',
			'text-align': 'center',
			'vertical-align': 'bottom'
		}),


		html.Div(
			children=[
			html.Label('Choose color map:'),
				dcc.Dropdown(
					id='colormap',
					options=[{'value': x, 'label': x} 
							 for x in colormaps],
					value='inferno',
					style={
						'width': '14.5vw'})
			],
			style={
			'display':'inline-block',
			'width': '14.5vw',
			'margin-left': '3vw', 
			'margin-right': '0vw',
			'margin-top': '2.5vh',
			'text-align': 'top',
			'vertical-align': 'bottom'
		}),


		html.Div(
			children=[
			html.Label('Choose method:'),
			dcc.Dropdown(
					id='method',
					options=[{'value': x, 'label': x} 
							 for x in methods],
					value=methods[0],
					style={
						'width': '17vw'})
			],
			style={
			'display': 'inline-block',
			'width': '17vw',
			'margin-right': '0vw',
			'margin-left':'2vw',
			'padding-left': '1vw',
			'margin-top': '2.5vh',
			'vertical-align': 'bottom'
		}),
			

		html.Div(
			children=[
			html.Label('Maximum iterations'),
				dcc.Input(
				id='steps',
				type='number',
				value=40,
				min=0,
				max=300,
				style={'margin-top': '0vh',
						'width': '15vw'})
			], 
			style={
			'display':'inline-block',
			'width': '15vw',
			'margin-left': '2vw', 
			'margin-right': '1vw',
			'margin-top': '2.5vh',
			'text-align': 'center',
			'vertical-align': 'top'
		}),

		
		html.Div(
			children=[
				html.Label('Specify Equation'),
					dcc.Input(
					id='equation',
					type='text',
					value='x^7.14-x-1',
					style={'margin-top': '1vh',
							'width': '20vw'})
			], 
			style={
			'display':'inline-block',
			'width': '20vw',
			'margin-left': '3vw', 
			'margin-right': '1vw',
			'margin-top': '2vh',
			'text-align': 'center',
			'vertical-align':'bottom'
		}),

		html.Button('Click to run', 
			id='button', 
			style={'display': 'inline-block',
					'margin-left': '2vw',
					'font-size': '1.4rem',
					'margin-top': '2vh',
					'vertical-align': 'bottom'}),
		html.Div(
			id='mathjax', 
			style={
				'textAlign': 'left',
				'font-family': "Open Sans", # "HelveticaNeue", "Helvetica Neue", Helvetica, Arial, sans-serif;', 
				'font-weight': 'normal',
				'margin-top': '2vh',
				'margin-left': '3vw',
				'font-size': '2.2rem',
				'display': 'inline-block',
				'margin-top': '4vh',
				'margin-bottom': '0vh',
				'vertical-align': 'bottom'
			}),


		html.Div(
			id='status', 
			style={
				'textAlign': 'left',
				'font-family': "Open Sans", # "HelveticaNeue", "Helvetica Neue", Helvetica, Arial, sans-serif;', 
				'font-weight': 'normal',
				'margin-top': '2vh',
				'margin-left': '3vw',
				'font-size': '2.2rem',
				'display': 'inline-block',
				'vertical-align': 'bottom'
			})
		],
		style={
			'display': 'inline-block',
			'width': '90vw',
			'margin-top': '0vh',
			'vertical-align': 'top'
		}
	),

	html.Img(
			id='image',
			style={'display': 'inline-block',
					'width': '90vw',
					'margin-left': '4vw',
					'margin-top': '3vh'}),

	dcc.Interval(id='trigger', interval=2000),
	# hidden div to store redis queue info
	html.Div(id='job', style={'display': 'none'}, children=dcc.Store(job)),

])

# responsive callbacks ('component_id' etc are not strictly necessary but 
# are included for clarity)
@app.callback(Output(component_id='image', component_property='src'),
			[Input(component_id='job', component_property='children'),
			Input(component_id='trigger', component_property='n_intervals'),
			Input(component_id='button', component_property='n_clicks')
			])
def update_redis(job, img, n_clicks):
	if n_clicks is None:
		raise PreventUpdate
		return 

	# wait while img is generated
	if q.count > 0:
		return 
	job_current = q.fetch_job('root_job')
	
	if not job_current:
		return

	# executes if redis job is complete
	if job_current.result:
		img = job_current.result
		return img


@app.callback(Output(component_id='job', component_property='children'),
			[Input(component_id='equation', component_property='value'),
			 Input(component_id='rbounds', component_property='value'),
			 Input(component_id='ibounds', component_property='value'),
			 Input(component_id='colormap', component_property='value'),
			 Input(component_id='steps', component_property='value'),
			 Input(component_id='method', component_property='value'),
			 Input(component_id='button', component_property='n_clicks')
			])
def display_juliaset(equation, rbounds, ibounds, colormap_value, steps_value, method, n_clicks):
	if n_clicks is None:
		raise PreventUpdate
	# do not update if redis queue has items waiting
	if q.count > 0:
		return ''
	# convert inputs to args and build array
	max_iterations = steps_value
	res_value = [1800, 1300]
	x_range = [i for i in rbounds.split(',')]
	y_range = [i for i in ibounds.split(',')]
	cmap = colormap_value

	if method == "Newton's":
		method = newton

	elif method == "Halley's":
		method = halley

	else:
		method = secant

	# send job to redis queue if not already there
	if not q.fetch_job('root_job'):
		q.enqueue(method, equation, max_iterations, x_range, y_range, res_value, cmap,
				ttl=1, failure_ttl=0.5, result_ttl=2, job_id='root_job')

	return ''

@app.callback(Output(component_id='button', component_property='n_clicks'),
		Input(component_id='image', component_property='src'))
def reset_clicks(img):
	if img and q.count == 0:
		return None
	else:
		return 1

@app.callback(
	Output(component_id='mathjax', component_property='children'),
	[Input(component_id='equation', component_property='value')])
def display_equation(equation):
	# convert equation to markdown and display
	final_string = ''
	i = 0
	while i in range(len(equation)):
		if equation[i] == '^':
			final_string += '^{'
			j = 0
			i += 1
			while equation[i+j] in '01234567890.':
				final_string += equation[i+j]
				j += 1
			final_string += '}'
			i = i + j
		else:
			final_string += equation[i]
			i += 1

	return f"\(  {final_string} \)"
	

@app.callback(
	Output(component_id='status', component_property='children'),
	[Input(component_id='button', component_property='n_clicks')])
def display_status(n_clicks):
	# show program status
	if n_clicks:
		return 'Running...'
	else:
		return ''

@app.callback(Output('trigger', 'interval'),
              [Input('image', 'src')])
def disable_interval(img):
    if img:
        return 60*60*1000 # one day
    else:
        return 2000


# run the app in the cloud
if __name__ == '__main__':
	# app.run_server(debug=True, port=8004)
	app.run_server(debug=True, host='0.0.0.0')




