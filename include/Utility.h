#pragma once

#include <cstdarg>

inline double min (int num, ...) {
	double temp, current_min;

	va_list arguments;
	va_start (arguments, num);
	current_min = va_arg (arguments, double);
	for (int x = 1; x < num; ++x) {
		temp = va_arg (arguments, double);
		if (temp < current_min) current_min = temp;
	}
	va_end (arguments);
	return current_min;
}
