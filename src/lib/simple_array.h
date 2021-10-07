/**
 * @file simple_array.h
 * @brief Functions to manipulate a simply array structure.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
*/

#ifndef _SIMPLEARRAY_
#define _SIMPLEARRAY_

#include <stdlib.h>
#include <string.h>

typedef struct sarray_ctx_
{
  void**  array;
  size_t* length;
  size_t  alength;
  size_t  esize;
  float   gfact;
} sarray_ctx;

inline static sarray_ctx* sarray_init(void** array, size_t *length, const size_t element_size, const float growing_factor)
{
  sarray_ctx* ctx=(sarray_ctx*)malloc(sizeof(sarray_ctx));
  *array=NULL;
  *length=0;
  ctx->array=array;
  ctx->length=length;
  ctx->alength=0;
  ctx->esize=element_size;
  ctx->gfact=growing_factor;
  return ctx;
}

inline static void sarray_clear( sarray_ctx* ctx)
{
  free(*ctx->array); *ctx->array=NULL;
  ctx->alength=0;
  *ctx->length=0;
}

inline static void sarray_free( sarray_ctx* ctx)
{
  sarray_clear(ctx);
  free(ctx);
}

inline static void sarray_clean_alloc(sarray_ctx* ctx)
{
  if(ctx->alength > *ctx->length) {
    ctx->alength=*ctx->length;
    *ctx->array=realloc(*ctx->array,ctx->alength*ctx->esize);
  }
}

inline static void sarray_empty(sarray_ctx* ctx)
{
  *ctx->length=0;
}

inline static int sarray_alloc(sarray_ctx* ctx, const size_t n)
{
  const size_t needed=*ctx->length+n;

  if(needed > ctx->alength) {
    ctx->alength=ctx->alength*ctx->gfact+1;

    if(needed > ctx->alength) ctx->alength=needed;
    *ctx->array=realloc(*ctx->array,ctx->alength*ctx->esize);

    if(!*ctx->array) return -1;
  }
  return 0;
}

inline static int sarray_grow(sarray_ctx* ctx, const size_t n)
{
  const size_t needed=*ctx->length+n;

  if(needed > ctx->alength) {
    ctx->alength=ctx->alength*ctx->gfact+1;

    if(needed > ctx->alength) ctx->alength=needed;
    *ctx->array=realloc(*ctx->array,ctx->alength*ctx->esize);

    if(!*ctx->array) return -1;
  }
  *ctx->length=needed;
  return 0;
}

inline static int sarray_grow_one(sarray_ctx* ctx)
{
  if(*ctx->length == ctx->alength) {
    ctx->alength=ctx->alength*ctx->gfact+1;
    *ctx->array=realloc(*ctx->array,ctx->alength*ctx->esize);

    if(!*ctx->array) return -1;
  }
  ++*ctx->length;
  return 0;
}

inline static void sarray_remove_last(sarray_ctx* ctx)
{
  if(*ctx->length > 0) --*ctx->length;
}

inline static int sarray_add_many_at(sarray_ctx* ctx, const size_t n, const size_t index)
{
  const size_t needed=*ctx->length+n;

  if(needed > ctx->alength) {
    ctx->alength=ctx->alength*ctx->gfact+1;

    if(needed > ctx->alength) ctx->alength=needed;
    *ctx->array=realloc(*ctx->array,ctx->alength*ctx->esize);

    if(!*ctx->array) return -1;
  }

  memmove(*ctx->array+index+n, *ctx->array+index, (*ctx->length-index)*ctx->esize);
  *ctx->length=needed;
  return 0;
}

inline static int sarray_add_space_at(sarray_ctx* ctx, const size_t index)
{
  if(*ctx->length == ctx->alength) {
    ctx->alength=ctx->alength*ctx->gfact+1;
    *ctx->array=realloc(*ctx->array,ctx->alength*ctx->esize);

    if(!*ctx->array) return -1;
  }
  memmove(*ctx->array+index+1, *ctx->array+index, (*ctx->length-index)*ctx->esize);
  ++*ctx->length;
  return 0;
}

inline static void sarray_remove_many_at(sarray_ctx* ctx, const size_t n, const size_t index)
{
  *ctx->length-=n;
  memmove(*ctx->array+index, *ctx->array+index+n, (*ctx->length-index)*ctx->esize);
}

inline static void sarray_remove_at(sarray_ctx* ctx, const size_t index)
{
  --*ctx->length;
  memmove(*ctx->array+index, *ctx->array+index+1, (*ctx->length-index)*ctx->esize);
}

#endif
