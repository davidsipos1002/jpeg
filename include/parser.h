#pragma once

#include <stdlib.h>
#include <string.h>

#include <format.h>
#include <misc.h>
#include <ds/utlist.h>

BEGIN_C_DECLS

// parse dht marker segment 
dht parse_dht(uint8_t *p);
void free_dht(dht t);


#define DUMP_QUANTIZATION
// extract quantization table
dqt *parse_dqt(uint8_t *p); 
void free_dqt(dqt *q);

END_C_DECLS
