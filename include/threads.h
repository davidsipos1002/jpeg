#pragma once

#include <pthread.h>
#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

typedef void *(thread_func) (void *);

pthread_t *thread_create(thread_func func, void *arg);
void thread_free(pthread_t *t);

pthread_mutex_t *mutex_create();
void mutex_free(pthread_mutex_t *m);

pthread_cond_t *cond_create();
void cond_free(pthread_cond_t *c);

END_C_DECLS
