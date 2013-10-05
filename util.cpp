//
//  util.cpp
//  AlgoUmMilotZuZerstoeren
//
//  Created by Martin Steinegger on 23.10.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include "util.h"
#include <iostream>
#include <stdlib.h>
#include <stdint.h>

void *memalign(size_t boundary, size_t size)
{
    void *pointer;
    if (posix_memalign(&pointer,boundary,size) != 0)
    {
        exit(3);
    }
    return pointer;
}



std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}


void strrev(char *p)
{
    char *q = p;
    while(q && *q) ++q;
    for(--q; p < q; ++p, --q)
        *p = *p ^ *q,
        *q = *p ^ *q,
        *p = *p ^ *q;
}