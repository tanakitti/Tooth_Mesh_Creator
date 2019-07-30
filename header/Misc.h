#pragma once

#include <string>
#include <vector>
#include <string.h>
#include <type_traits>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include "buildparams.h"

namespace pd
{

    template <typename T>
    static
    T clamp(const T& n, const T& lower, const T& upper) {
#if __cplusplus >= 201703L
        return std::clamp<T>( n, lower, upper );
#else
        return std::max(lower, std::min(n, upper));
#endif //__cplusplus
    }

    template <typename T>
    static
    T lerp(const T &t, const T &low, const T &high)
    {
        return low + t*(high - low);
    }

    template <typename T>
    static
    T unlerp(const T &t, const T &low, const T &high)
    {
#if _DEBUG
        T denom = high - low;
	if (denom == 0)
		return 0;
#endif
        return (t - low) / (high - low);
    }

    template <typename T>
    static
    T remap(const T &t, const T &low1, const T &high1, const T &low2, const T &high2)
    {
        return lerp<T>(unlerp<T>(t, low1, high1), low2, high2);
    }

    static
    std::vector<std::string> Split( std::string &str, const char *delims )
    {
        std::vector<std::string> result;
        char *dup = strdup( str.c_str() );
        char *tok = strtok( dup, delims );

        while (tok != nullptr) {
            result.push_back( std::string(tok) );
            tok = strtok( nullptr, delims );
        }

        return result;
    }

    static
    std::string GetVersionDate()
    {
        return std::string( __DATE__ );
    }

    static
    std::string GetVersionTime()
    {
        return std::string( __TIME__ );
    }

    static
    std::string GetVersionDateTime()
    {
        return GetVersionDate() + " " + GetVersionTime();
    }

    static
    std::string GetVersionGitCommit()
    {
        return std::string( BUILD_GIT_COMMIT );
    }

    static
    std::string GetVersionGitBranch()
    {
        return std::string( BUILD_GIT_BRANCH );
    }

    static
    std::string GetRuntimeDate()
    {
        static std::string RuntimeDate;
        static bool once = [](){
            auto t = std::time(nullptr);
            auto tm = *std::localtime(&t);
            std::ostringstream date;
            date << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
            RuntimeDate = date.str();
            return true;
        } ();
        return RuntimeDate;
    }

    static
    void WriteVersionInfo( std::ofstream &file )
    {
        file << "# Runtime " << GetRuntimeDate() << std::endl;
        file << "# Compiletime " << GetVersionDate() << std::endl;
        file << "# Gitcommit " << GetVersionGitCommit() << std::endl;
        file << "# Gitbranch " << GetVersionGitBranch() << std::endl;
    }

}