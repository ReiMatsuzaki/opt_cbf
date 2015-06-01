#ifndef MACROS_HPP
#define MACROS_HPP

#include <stdio.h>

#define SUB_LOCATION(msg) msg=__FILE__;\
                          msg+=":";\
                          msg+=__FUNCTION__;\
                          msg+=":";\
			  char line[10];\
			  sprintf(line, "%d", __LINE__);\
                          msg+=line;
#endif


