#ifndef _INFINDIVIDUAL_
#define _INFINDIVIDUAL_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

struct infindividual
{
    struct infindividual* next;
    void* dataptr;
    float* etimes;
    float* etimesptr;
    int nevents;
};

struct iillist
{
  struct infindividual* head;
  struct infindividual* tail;
  struct infindividual* recycling_head;
};

static inline void iillist_init(struct iillist* list){memset(list, 0, sizeof(struct iillist));}

static inline void iillist_add_element(struct iillist* list, struct infindividual* element){element->next=list->head; list->head=element; if(!list->tail) list->tail=element;}
static inline void iillist_add_many_elements(struct iillist* list, struct infindividual* head, struct infindividual* tail){tail->next=list->head; list->head=head; list->tail=tail;}

static inline void iillist_recycle_element(struct iillist* list, struct infindividual* prev_element, struct infindividual* element){if(prev_element) prev_element->next=element->next; else list->head=element->next; element->next=list->recycling_head; list->recycling_head=element; if(list->tail==element) list->tail=prev_element;}
static inline void iillist_recycle_detached_element(struct iillist* list, struct infindividual* element){element->next=list->recycling_head; list->recycling_head=element;}

static inline void iillist_recycle_all(struct iillist* list){list->tail->next=list->recycling_head; list->recycling_head=list->head; list->head=list->tail=NULL;}
static inline void iillist_clear_recycling(struct iillist* list){struct infindividual* cur; while(list->recycling_head) {cur=list->recycling_head; list->recycling_head=cur->next; free(cur);}}
static inline void iillist_clear(struct iillist* list){struct infindividual* cur; while(list->head) {cur=list->head; list->head=cur->next; free(cur);} iillist_clear_recycling(list);}

static inline struct infindividual* iillist_new_element(struct iillist* list){struct infindividual* ret; if(list->recycling_head) {ret=list->recycling_head; list->recycling_head=ret->next; return ret;} ret=(struct infindividual*)malloc(sizeof(struct infindividual)); memset(ret, 0, sizeof(struct infindividual)); return ret;} 

#endif
