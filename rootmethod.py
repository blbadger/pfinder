# rootmethod.py

# import standard libraries
import io
import base64

# impot third party libraries
import numpy as np 
import matplotlib.pyplot as plt 
from CalculateFaster import OptiCalculate
from CoCalculate import ComplexCalculate
import numexpr as ne


def convert_to_binary(arr, cmap):
	buf = io.BytesIO()
	plt.style.use('dark_background')
	plt.imsave(buf, arr, cmap=cmap, format='png')
	data = base64.b64encode(buf.getbuffer()).decode("utf8") 
	return "data:image/png;base64,{}".format(data)


def newton(equation, max_iterations, x_range, y_range, res_value, cmap):
	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	# initialization of Calculate object
	nondiff = ComplexCalculate(equation, differentiate=False)
	diff = ComplexCalculate(equation, differentiate=True)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z)
		z_array = z_array - f_now / f_prime_now

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		
	arr = iterations_until_rooted
	return convert_to_binary(arr, cmap)



def newton_optimized(equation, max_iterations, x_range, y_range, res_value, cmap):
	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	# initialization of OptiCalculate object
	nondiff = OptiCalculate(equation, differentiate=False)
	diff = OptiCalculate(equation, differentiate=True)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z)
		z_array = ne.evaluate('z_array - f_now / f_prime_now')

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		
	arr = iterations_until_rooted
	return convert_to_binary(arr, cmap)



def halley(equation, max_iterations, x_range, y_range, res_value, cmap):
	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	nondiff = ComplexCalculate(equation, differentiate=False)
	diff = ComplexCalculate(equation, differentiate=True)

	# double derivative initialization
	diff_stirng = diff.to_string()
	double_diff = ComplexCalculate(diff_string, differentiate=True)


	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z) # first derivative evaluation
		f_double_prime_now = double_diff.evaluate(z) # second derivative evaluation
		z_array = z - (2*f_now * f_prime_now / (2*(f_prime_now)**2 - f_now * f_double_prime_now))

		# test the boolean map for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.000000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	arr = iterations_until_rooted
	return convert_to_binary(arr, cmap)



def halley_optimized(equation, max_iterations, x_range, y_range, res_value, cmap):
	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	#
	nondiff = OptiCalculate(equation, differentiate=False)
	diff = OptiCalculate(equation, differentiate=True)

	# double derivative initialization
	diff_stirng = diff.to_string()
	double_diff = OptiCalculate(diff_string, differentiate=True)


	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z) # first derivative evaluation
		f_double_prime_now = double_diff.evaluate(z) # second derivative evaluation
		z_array = ne.evaluate('z - (2*f_now * f_prime_now / (2*(f_prime_now)**2 - f_now * f_double_prime_now))')

		# test the boolean map for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.000000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	arr = iterations_until_rooted
	return convert_to_binary(arr, cmap)



def secant(equation, max_iterations, x_range, y_range, res_value, cmap):
	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j
	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000
	zeros = np.zeros(z_array.shape) 
	z_0 = (z_array - zeros)/2 # setting the initial guess to half the distance to the origin from the second guess, which is plotted

	# initializatoin of ComplexCalculate object
	nondiff = ComplexCalciulate(equation, differentiate=False)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_previous = nondiff.evaluate(z_0)
		f_now = nondiff.evaluate(z)
		z_array = z - f_now * (z - z_0)/(f_now - f_previous)

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		z_0 = z 


	arr = iterations_until_rooted
	return convert_to_binary(arr, cmap)


def secant_optimized(equation, max_iterations, x_range, y_range, res_value, cmap):
	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j
	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000
	zeros = np.zeros(z_array.shape) 
	z_0 = ne.evaluate('(z_array - zeros)/2') # setting the initial guess to half the distance to the origin from the second guess, which is plotted

	# initialize OptiCalculate object
	nondiff = OptiCalculate(equation, differentiate=False)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_previous = nondiff.evaluate(z_0)
		f_now = nondiff.evaluate(z)
		z_array = ne.evaluate('z - f_now * (z - z_0)/(f_now - f_previous)')

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		z_0 = z 


	arr = iterations_until_rooted
	return convert_to_binary(arr, cmap)