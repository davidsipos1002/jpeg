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

// extract quantization table
dqt *parse_dqt(uint8_t *p, uint16_t *l); 
void free_dqt(dqt *q);

// parse dnl marker segment
uint16_t parse_dnl(uint8_t *p, uint16_t *l); 

// parse dri marker segment
uint16_t parse_dri(uint8_t *p, uint16_t *l);

// parse frame header
fhdr *parse_fhdr(uint8_t *p, uint16_t *l);
void free_fhdr(fhdr *fhdr);

// parse scan header
shdr *parse_shdr(uint8_t *p, uint16_t *l);
void free_shdr(shdr *s);

// verify JFIF
uint8_t parse_app0(uint8_t *p, uint16_t *l);

uint16_t get_marker_seg_len(uint8_t *p);

END_C_DECLS
