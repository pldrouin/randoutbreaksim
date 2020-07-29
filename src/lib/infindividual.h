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
  struct infindividual* rec_head;
};

static inline void iillist_init(struct iillist* list){memset(list, 0, sizeof(struct iillist));}

static inline void iillist_addelement(struct iillist* list, struct infindividual* element){element->next=list->head; list->head=element; if(!list->tail) list->tail=element;}
static inline void iillist_addmanyelements(struct iillist* list, struct infindividual* head, struct infindividual* tail){tail->next=list->head; list->head=head; list->tail=tail;}

static inline void iillist_recelement(struct iillist* list, struct infindividual* prev_element, struct infindividual* element){if(prev_element) prev_element->next=element->next; else list->head=element->next; element->next=list->rec_head; list->rec_head=element; if(list->tail==element) list->tail=prev_element;}

static inline void iillist_rec_all(struct iillist* list){list->tail->next=list->rec_head; list->rec_head=list->head; list->head=list->tail=NULL;}
static inline void iillist_clear_rec(struct iillist* list){struct infindividual* cur; while(list->rec_head) {cur=list->rec_head; list->rec_head=cur->next; free(cur);}}
static inline void iillist_clear(struct iillist* list){struct infindividual* cur; while(list->head) {cur=list->head; list->head=cur->next; free(cur);} iillist_clear_rec(list);}

static inline struct infindividual* iillist_newelement(struct iillist* list){struct infindividual* ret; if(list->rec_head) {ret=list->rec_head; list->rec_head=ret->next; return ret;} ret=(struct infindividual*)malloc(sizeof(struct infindividual)); memset(ret, 0, sizeof(struct infindividual)); return ret;} 

#endif
