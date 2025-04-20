#ifndef MY_RAND_H
#define MY_RAND_H

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

static int my_rand(int range)
{

    return random()%range;
}

#endif
