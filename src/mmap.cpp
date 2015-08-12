#include <filesystem/path.h>
#include <layer/mmap.h>
#include <layer/log.h>

#if defined(__LINUX__) || defined(__OSX__)
# include <sys/mman.h>
# include <fcntl.h>
#elif defined(__WINDOWS__)
# include <windows.h>
#endif

NAMESPACE_BEGIN(layer)

struct MemoryMappedFile::MemoryMappedFilePrivate {
	fs::path filename;
#if defined(__WINDOWS__)
	HANDLE file;
	HANDLE fileMapping;
#endif
	size_t size;
	void *data;
	bool readOnly;
	bool temp;

	MemoryMappedFilePrivate(const fs::path &f = "", size_t s = 0)
		: filename(f), size(s), data(NULL), readOnly(false), temp(false) {}

	void create() {
		#if defined(__LINUX__) || defined(__OSX__)
			int fd = open(filename.str().c_str(), O_RDWR | O_CREAT | O_TRUNC, 0664);
			if (fd == -1)
				Error("Could not open \"%s\"!", filename.str());
			int result = lseek(fd, size-1, SEEK_SET);
			if (result == -1)
				Error("Could not set file size of \"%s\"!", filename.str());
			result = write(fd, "", 1);
			if (result != 1)
				Error("Could not write to \"%s\"!", filename.str());
			data = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
			if (data == NULL)
				Error("Could not map \"%s\" to memory!", filename.str());
			if (close(fd) != 0)
				Error("close(): unable to close file!");
		#elif defined(__WINDOWS__)
			file = CreateFileW(filename.wstr().c_str(), GENERIC_WRITE | GENERIC_READ,
				FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
			if (file == INVALID_HANDLE_VALUE)
				Error("Could not open \"%s\": %s", filename.str(),
					lastErrorText());
			fileMapping = CreateFileMapping(file, NULL, PAGE_READWRITE, 0,
				static_cast<DWORD>(size), NULL);
			if (fileMapping == NULL)
				Error("CreateFileMapping: Could not map \"%s\" to memory: %s",
					filename.str(), lastErrorText());
			data = (void *) MapViewOfFile(fileMapping, FILE_MAP_WRITE, 0, 0, 0);
			if (data == NULL)
				Error("MapViewOfFile: Could not map \"%s\" to memory: %s",
					filename.str(), lastErrorText());
		#endif
		readOnly = false;
	}

	void createTemp() {
		readOnly = false;
		temp = true;

		#if defined(__LINUX__) || defined(__OSX__)
			char *path = strdup("/tmp/mmap_XXXXXX");
			int fd = mkstemp(path);
			if (fd == -1)
				Error("Unable to create temporary file (1): %s", strerror(errno));
			filename = path;
			free(path);

			int result = lseek(fd, size-1, SEEK_SET);
			if (result == -1)
				Error("Could not set file size of \"%s\"!", filename.str());
			result = write(fd, "", 1);
			if (result != 1)
				Error("Could not write to \"%s\"!", filename.str());

			data = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
			if (data == NULL)
				Error("Could not map \"%s\" to memory!", filename.str());

			if (close(fd) != 0)
				Error("close(): unable to close file!");
		#elif defined(__WINDOWS__)
			WCHAR tempPath[MAX_PATH];
			WCHAR tempFilename[MAX_PATH];

			unsigned int ret = GetTempPathW(MAX_PATH, tempPath);
			if (ret == 0 || ret > MAX_PATH)
				SError("GetTempPath failed(): %s", lastErrorText().c_str());

			ret = GetTempFileNameW(tempPath, L"mitsuba", 0, tempFilename);
			if (ret == 0)
				SError("GetTempFileName failed(): %s", lastErrorText().c_str());

			file = CreateFileW(tempFilename, GENERIC_READ | GENERIC_WRITE,
				0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);

			if (file == INVALID_HANDLE_VALUE)
				Error("Error while trying to create temporary file: %s",
					lastErrorText().c_str());

			filename = fs::path(tempFilename);

			fileMapping = CreateFileMapping(file, NULL, PAGE_READWRITE, 0,
				static_cast<DWORD>(size), NULL);
			if (fileMapping == NULL)
				Error("CreateFileMapping: Could not map \"%s\" to memory: %s",
					filename.str(), lastErrorText());
			data = (void *) MapViewOfFile(fileMapping, FILE_MAP_WRITE, 0, 0, 0);
			if (data == NULL)
				Error("MapViewOfFile: Could not map \"%s\" to memory: %s",
					filename.str(), lastErrorText());
		#endif
	}

	void map() {
		if (!filename.exists())
			Error("The file \"%s\" does not exist!", filename.str());
		size = filename.file_size();

		#if defined(__LINUX__) || defined(__OSX__)
			int fd = open(filename.str().c_str(), readOnly ? O_RDONLY : O_RDWR);
			if (fd == -1)
				Error("Could not open \"%s\"!", filename.str());
			data = mmap(NULL, size, PROT_READ | (readOnly ? 0 : PROT_WRITE), MAP_SHARED, fd, 0);
			if (data == NULL)
				Error("Could not map \"%s\" to memory!", filename.str());
			if (close(fd) != 0)
				Error("close(): unable to close file!");
		#elif defined(__WINDOWS__)
			file = CreateFileW(filename.wstr().c_str(), GENERIC_READ | (readOnly ? 0 : GENERIC_WRITE),
				FILE_SHARE_WRITE|FILE_SHARE_READ, NULL, OPEN_EXISTING,
				FILE_ATTRIBUTE_NORMAL, NULL);
			if (file == INVALID_HANDLE_VALUE)
				Error("Could not open \"%s\": %s", filename.str(),
					lastErrorText());
			fileMapping = CreateFileMapping(file, NULL, readOnly ? PAGE_READONLY : PAGE_READWRITE, 0, 0, NULL);
			if (fileMapping == NULL)
				Error("CreateFileMapping: Could not map \"%s\" to memory: %s",
					filename.str(), lastErrorText());
			data = (void *) MapViewOfFile(fileMapping, readOnly ? FILE_MAP_READ : FILE_MAP_WRITE, 0, 0, 0);
			if (data == NULL)
				Error("MapViewOfFile: Could not map \"%s\" to memory: %s",
					filename.str(), lastErrorText());
		#endif
	}

	void unmap() {
		Trace("Unmapping \"%s\" from memory", filename.str());

		#if defined(__LINUX__) || defined(__OSX__)
			if (temp) {
				/* Temporary file that will be deleted in any case:
				   invalidate dirty pages to avoid a costly flush to disk */
				int retval = msync(data, size, MS_INVALIDATE);
				if (retval != 0)
					Error("munmap(): unable to unmap memory: %s", strerror(errno));
			}

			int retval = munmap(data, size);
			if (retval != 0)
				Error("munmap(): unable to unmap memory: %s", strerror(errno));
		#elif defined(__WINDOWS__)
			if (!UnmapViewOfFile(data))
				Error("UnmapViewOfFile(): unable to unmap memory: %s", lastErrorText().c_str());
			if (!CloseHandle(fileMapping))
				Error("CloseHandle(): unable to close file mapping: %s", lastErrorText().c_str());
			if (!CloseHandle(file))
				Error("CloseHandle(): unable to close file: %s", lastErrorText().c_str());
		#endif

		if (temp) {
			try {
				filename.remove_file();
			} catch (...) {
				Warn("unmap(): Unable to delete file \"%s\"", filename.str());
			}
		}

		data = NULL;
		size = 0;
	}
};

MemoryMappedFile::MemoryMappedFile()
	: d(new MemoryMappedFilePrivate()) { }

MemoryMappedFile::MemoryMappedFile(const fs::path &filename, size_t size)
	: d(new MemoryMappedFilePrivate(filename, size)) {
	Trace("Creating memory-mapped file \"%s\" (%s)..",
		filename.str(), memString(d->size).c_str());
	d->create();
}

MemoryMappedFile::MemoryMappedFile(const fs::path &filename, bool readOnly)
	: d(new MemoryMappedFilePrivate(filename)) {
	d->readOnly = readOnly;
	d->map();
    Trace("Mapped \"%s\" into memory (%s)..",
          filename.str(), memString(d->size).c_str());
}

MemoryMappedFile::~MemoryMappedFile() {
	if (d->data) {
		try {
			d->unmap();
		} catch (std::exception &e) {
			/* Don't throw exceptions from a destructor */
			Warn("%s", e.what());
		}
	}
}

void MemoryMappedFile::resize(size_t size) {
	if (!d->data)
		Error("Internal error in MemoryMappedFile::resize()!");
	bool temp = d->temp;
	d->temp = false;
	d->unmap();
	d->filename.resize_file(size);
	d->size = size;
	d->map();
	d->temp = temp;
}

void *MemoryMappedFile::data() {
	return d->data;
}

/// Return a pointer to the file contents in memory (const version)
const void *MemoryMappedFile::data() const {
	return d->data;
}

size_t MemoryMappedFile::size() const {
	return d->size;
}

bool MemoryMappedFile::readOnly() const {
	return d->readOnly;
}

const fs::path &MemoryMappedFile::filename() const {
	return d->filename;
}

MemoryMappedFile *MemoryMappedFile::createTemporary(size_t size) {
	MemoryMappedFile *result = new MemoryMappedFile();
	result->d->size = size;
	result->d->createTemp();
	return result;
}

std::string MemoryMappedFile::toString() const {
	std::ostringstream oss;
	oss << "MemoryMappedFile[filename=\""
		<< d->filename.str() << "\", size="
		<< memString(d->size) << "]";
	return oss.str();
}

NAMESPACE_END(layer)
