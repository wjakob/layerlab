/*
    mmap.h -- Portable memory mapped file implementation

    Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <layer/common.h>
#include <memory>

NAMESPACE_BEGIN(layer)

/// Basic cross-platform abstraction for memory mapped files
class MemoryMappedFile {
public:
	/// Create a new memory-mapped file of the specified size
	MemoryMappedFile(const fs::path &filename, size_t size);

	/// Map the specified file into memory
	MemoryMappedFile(const fs::path &filename, bool readOnly = true);

	/// Return a pointer to the memory-mapped file contents
	void *data();

	/// Return a pointer to the memory-mapped file contents (const version)
	const void *data() const;

	/// Return the size of the mapped region
	size_t size() const;

	/**
	 * \brief Resize the memory-mapped file
	 *
	 * This involves remapping the file, which will
	 * generally change the pointer obtained via data()
	 */
	void resize(size_t size);

	/// Return the associated filename
	const fs::path &filename() const;

	/// Return whether the mapped memory region is read-only
	bool readOnly() const;

	/// Return a string representation
	std::string toString() const;

	/**
	 * \brief Create a temporary memory-mapped file
	 *
	 * \remark When closing the mapping, the file is
	 * automatically deleted.
	 */
	static MemoryMappedFile* createTemporary(size_t size);

	/// Release all resources
	virtual ~MemoryMappedFile();

protected:
	/// Internal constructor
	MemoryMappedFile();

private:
	struct MemoryMappedFilePrivate;
	std::unique_ptr<MemoryMappedFilePrivate> d;
};

NAMESPACE_END(layer)
