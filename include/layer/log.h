/*
    common.h -- Basic macros and type definitions

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/common.h>

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4702) // warning C4702: unreachable code
#endif
#include <tinyformat.h>
#if defined(_MSC_VER)
#pragma warning(pop)
#endif


#define Error(...) \
    throw std::runtime_error(tfm::format(__VA_ARGS__))

#define Trace(format, ...) \
    tfm::printf("Trace: " format "\n", ##__VA_ARGS__)

#define Warn(format, ...) \
    tfm::printf("Warning: " format "\n", ##__VA_ARGS__)

#define Log(format, ...) \
    tfm::printf("Log: " format "\n", ##__VA_ARGS__)

NAMESPACE_BEGIN(layer)

/**
 * \brief Convert a time difference (in seconds) into a human-readable string
 *        representation
 * \param time
 *        Time difference in (fractional) sections
 * \param precise
 *        When set to true, a higher-precision string representation is
 *        generated.
 */
extern std::string timeString(Float time, bool precise = false);

/**
 * \brief Convert a memory amount (in bytes) into a human-readable string
 *        representation
 * \param size
 *        An unsigned integer size value
 * \param precise
 *        When set to true, a higher-precision string representation is
 *        generated.
 */
extern std::string memString(size_t size, bool precise = false);

#if defined(__WINDOWS__)
/// Return a string version of GetLastError()
extern std::string lastErrorText();
#endif

NAMESPACE_END(layer)
