#pragma once

const size_t BYTES_IN_KILOBYTE = 1024;
const size_t BYTES_IN_MEGABYTE = 1024*1024;
const size_t BYTES_IN_GIGABYTE = 1024*1024*1024;


size_t bytesToKb(size_t bytes){
	return bytes/1024;
}

size_t bytesToMb(size_t bytes){
	return bytes/1024/1024;
}

size_t bytesToGb(size_t bytes){
	return bytes/1024/1024/1024;
}

std::string memoryString(size_t sizeInBytes){
	std::stringstream s;
	if(sizeInBytes<BYTES_IN_KILOBYTE) s << sizeInBytes << " bytes";
	else if(sizeInBytes<BYTES_IN_MEGABYTE){
		s << bytesToKb(sizeInBytes) << "Kb";
	}
	else if(sizeInBytes<BYTES_IN_GIGABYTE) {
		size_t residue = sizeInBytes % BYTES_IN_MEGABYTE;
		s << bytesToMb(sizeInBytes) << "Mb  " << bytesToKb(residue) << "Kb";
	}
	else {
		size_t residue = sizeInBytes % BYTES_IN_GIGABYTE;//TODO make shift
		s << bytesToGb(sizeInBytes) << "Gb  " << bytesToMb(residue) << "Mb";
	}
	return s.str();
}

size_t TOTAL_MEMORY_ALLOCATED_IN_BYTES = 0;

void trackMemory(size_t sizeInBytes){
	TOTAL_MEMORY_ALLOCATED_IN_BYTES+=sizeInBytes;
    std::cout << "memory allocated: " << memoryString(sizeInBytes) << std::endl;
	std::cout << "total memory allocated: " << memoryString(TOTAL_MEMORY_ALLOCATED_IN_BYTES) << std::endl;

}



