#ifndef _JUNCTION_H_
#define _JUNCTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>
#include "utils.h"
#include "uthash.h"

static inline int min_mismatch(char* str, char* pattern){
	if(str == NULL || pattern == NULL) die("[%s] input error"); 
	register int i, j, n;
	int min_mis_match = strlen(pattern)+1;
	char substring[strlen(pattern)+1];
	for(i=0;i<strlen(str)-strlen(pattern);i++){
		n = 0;
		memset(substring, '\0', sizeof(substring));
		strncpy(substring, &str[i], strlen(pattern));
		for(j=0; j<strlen(pattern); j++){if(toupper(pattern[j]) != toupper(substring[j])){n++;}}
		if (n < min_mis_match){min_mis_match = n;} // update min_mis_match
	}
	return min_mis_match;
} 

#endif
