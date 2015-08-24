#include <layer/log.h>
#include <cmath>
#if defined(__WINDOWS__)
#include <windows.h>
#endif

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

#if defined(__WINDOWS__)
std::string lastErrorText() {
    DWORD errCode = GetLastError();
    char *errorText = NULL;
    if (!FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER
        | FORMAT_MESSAGE_FROM_SYSTEM
        | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        errCode,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR)&errorText,
        0,
        NULL)) {
        return "Internal error while looking up an error code";
    }
    std::string result(errorText);
    LocalFree(errorText);
    return result;
}
#endif

NAMESPACE_END(layer)
