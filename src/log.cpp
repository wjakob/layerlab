#include <layer/log.h>
#include <cmath>

NAMESPACE_BEGIN(layer)

std::string timeString(Float time, bool precise) {
	if (std::isnan(time))
		return "nan";
	else if (std::isinf(time))
		return "inf";

	char suffix = 's';
	if (time > 60) {
		time /= 60; suffix = 'm';
		if (time > 60) {
			time /= 60; suffix = 'h';
			if (time > 12) {
				time /= 12; suffix = 'd';
			}
		}
	}

	std::ostringstream os;
	os.precision(precise ? 4 : 1);
	os << std::fixed << time << suffix;

	return os.str();
}

std::string memString(size_t size, bool precise) {
	Float value = (Float) size;
	const char *suffixes[] = {
		"B", "KiB", "MiB", "GiB", "TiB", "PiB"
	};
	int suffix = 0;
	while (suffix < 5 && value > 1024.0f) {
		value /= 1024.0f; ++suffix;
	}

	std::ostringstream os;
	os.precision(suffix == 0 ? 0 : (precise ? 4 : 1));
	os << std::fixed << value << " " << suffixes[suffix];

	return os.str();
}

NAMESPACE_END(layer)
