#include <threads.h>


pthread_t *thread_create(thread_func func, void *arg)
{
    pthread_t *t;
    safeMalloc(t, sizeof(pthread_t));
    pthread_create(t, NULL, func, arg);
    return t;
}

void thread_free(pthread_t *t)
{
    pthread_join(*t, NULL);
    free(t);
}

pthread_mutex_t *mutex_create()
{
    pthread_mutex_t *m;
    safeMalloc(m, sizeof(pthread_mutex_t));
    pthread_mutex_init(m, NULL);
    return m;
}

void mutex_free(pthread_mutex_t *m)
{
    pthread_mutex_destroy(m);
    free(m);
}

pthread_cond_t *cond_create()
{
    pthread_cond_t *c;
    safeMalloc(c, sizeof(pthread_cond_t));
    pthread_cond_init(c, NULL);
    return c;
}

void cond_free(pthread_cond_t *c)
{
    pthread_cond_destroy(c);
    free(c);
}
