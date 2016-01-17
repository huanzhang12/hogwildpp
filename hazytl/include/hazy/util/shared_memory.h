// Copyright 2012 Victor Bittorf, Chris Re
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Hazy Template Library 
// Created: October 2012
// Primary Author: Victor Bittorf (bittorf [at] cs.wisc.edu)


#ifndef HAZY_UTIL_SHARED_MEMORY_H
#define HAZY_UTIL_SHARED_MEMORY_H

#include <sys/types.h>
#include <sys/ipc.h> 
#include <sys/shm.h>

namespace hazy {
namespace util {

/*! \brief Gets the memory segment and returns a pointer
 * The pointer must be detached after use.
 * \note this can fail for many reasons, will print out error codes.
 * \param key shared memory segment's key
 * \param size the size in bytes of the segment
 * \return the pointer to the memory block or NULL if it failed.
 */
char* AttachSharedSegment(key_t key, size_t size) {
  int shmflag; // permission flag
  int shmid; // memory segment id
  char* shmaddr; //pointer to the head of the segment

  shmflag = IPC_CREAT | 0666;

  if ((shmid = shmget(key, size, shmflag)) < 0) {
    // [EACCES]
    // A shared memory identifier exists for key but operation permission as
    // specified by the low-order nine bits of shmflg would not be granted. See
    // IPC.
    // [EEXIST]
    // A shared memory identifier exists for the argument key but
    // (shmflg&IPC_CREAT)&&(shmflg&IPC_EXCL) is non-zero.
    // [EINVAL]
    // The value of size is less than the system-imposed minimum or greater
    // than the system-imposed maximum, or a shared memory identifier exists
    // for the argument key but the size of the segment associated with it is
    // less than size and size is not 0.
    // [ENOENT]
    // A shared memory identifier does not exist for the argument key and
    // (shmflg&IPC_CREAT) is 0.
    // [ENOMEM]
    // A shared memory identifier and associated shared memory segment are to
    // be created but the amount of available physical memory is not sufficient
    // to fill the request.
    // [ENOSPC]
    // A shared memory identifier is to be created but the system-imposed limit
    // on the maximum number of allowed shared memory identifiers system-wide
    // would be exceeded.
    perror("shmget failed");
    return NULL;
  }

  // use the shared identifier to attached the segment to our address space
  shmaddr = static_cast<char*>(shmat(shmid, NULL, 0));
  if(shmaddr == (char *)-1) {
    // [EACCES]
    // Operation permission is denied to the calling process, see IPC.
    // [EINVAL]
    // The value of shmid is not a valid shared memory identifier; the shmaddr
    // is not a null pointer and the value of
    // (shmaddr-((ptrdiff_t)shmaddr%SHMLBA)) is an illegal address for
    // attaching shared memory; or the shmaddr is not a null pointer,
    // (shmflg&SHM_RND) is 0 and the value of shmaddr is an illegal address for
    // attaching shared memory.
    // [EMFILE]
    // The number of shared memory segments attached to the calling process
    // would exceed the system-imposed limit.
    // [ENOMEM]
    // The available data space is not large enough to accommodate the shared
    // memory segment.
    perror("shmat failed");
    return NULL;
  }
  return shmaddr;
}

/*! \brief Detaches from the shared segment given by the attached address
 * \param shmaddr an addressed returned by AttachSharedSegment(...)
 * \return true if successful, on error returns false
 */
bool DetachSharedSegment(char* shmaddr) {
  if (shmdt(shmaddr) < 0) {
    // Only 1 error condition
    // [EINVAL]
    // The value of shmaddr is not the data segment start address of a shared
    // memory segment.
    perror("shmdat failed");
    return false;
  }
  return true;
}

} //namespace util
} //namespace hazy
#endif
