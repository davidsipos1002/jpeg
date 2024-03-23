#pragma once

#include <stdlib.h>

#include <format.h>
#include <misc.h>
#include <ds/utlist.h>

BEGIN_C_DECLS

// parse dht marker segment 
dht parse_dht(uint8_t *p);
void free_dht(dht t);

END_C_DECLS
