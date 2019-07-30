#pragma once

#ifndef _WIN32
#include <mutex>
#else
#include "Engine.h"
#endif


/**
 * Wrapper class for mutex and critical section.
 * On Windows, it uses a critical section.
 * On other systems, mutex is used.
 *
 * @author Kevin Marnholz
 */
class FMutex
{

public:

    void Lock()
    {
#ifdef _WIN32
        CriticalSection.Lock();
#else
        Mutex.lock();
#endif
    }

    bool TryLock()
    {
#ifdef _WIN32
        return CriticalSection.TryLock();
#else
        return Mutex.try_lock();
#endif
    }

    void Unlock()
    {
#ifdef _WIN32
        CriticalSection.Unlock();
#else
        Mutex.unlock();
#endif
    }




private:


#ifdef _WIN32
    FCriticalSection CriticalSection;
#else
    std::mutex Mutex;
#endif

};