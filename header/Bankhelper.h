//
// Created by Tanakit Sachati on 2019-07-15.
//

#ifndef UNTITLED4_BANKHELPER_H
#define UNTITLED4_BANKHELPER_H

#include <iostream>
#include <sstream>

template <typename T>
void print(T t)
{
    std::cout << t << std::endl;
}

template<typename T, typename... Args>
void print(T t, Args... args)
{
    std::cout << t << " ";
    print(args...) ;
}

template <typename T>
void printArray(T arr,int size,float great){
    for ( int i = 0; i < size; i++ )
        if(arr [ i ]>great){
            std::cout << i<< " : "<< arr [ i ] << std::endl;
        }

}

template <typename T>
std::string toString (T arg)
{
    std::stringstream ss;
    ss << arg;
    return ss.str ();
}

#endif //UNTITLED4_BANKHELPER_H
