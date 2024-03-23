#pragma once

#include <stdlib.h>
#include <string.h>

#include <jpeg/format.h>
#include <misc.h>
#include <ds/utlist.h>

BEGIN_C_DECLS

// parse dht marker segment 
dht parse_dht(uint8_t *p, uint16_t *l);
void free_dht(dht t);

#define DUMP_QUANTIZATION
// extract quantization table
dqt *parse_dqt(uint8_t *p, uint16_t *l); 
void free_dqt(dqt *q);

// parse dnl marker segment
uint16_t parse_dnl(uint8_t *p, uint16_t *l); 

// parse dri marker segment
uint16_t parse_dri(uint8_t *p, uint16_t *l);

#define DUMP_FRAME_HEADER
// parse frame header
fhdr *parse_fhdr(uint8_t *p, uint16_t *l);
void free_fhdr(fhdr *fhdr);

#define DUMP_SCAN_HEADER
// parse scan header
shdr *parse_shdr(uint8_t *p, uint16_t *l);
void free_shdr(shdr *s);

uint16_t get_marker_seg_len(uint8_t *p);

END_C_DECLS
