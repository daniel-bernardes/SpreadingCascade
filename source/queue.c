#include <stdlib.h>
#include <assert.h>

typedef struct _Queue {
  int size;
  int begin;
  int end;
  int *nodes;
} Queue;

inline int queue_empty(Queue *q){ return (q->begin == q->end); }
inline int queue_full(Queue *q) { return (q->begin == (q->end+1) % q->size); }

Queue *queue_new(int size) {
  Queue *q = (Queue *) malloc(sizeof(Queue));
  assert(q != NULL);
  q->size = size+1;
  q->nodes = (int *) malloc(q->size * sizeof(int));
  assert(q->nodes != NULL);
  q->begin = q->end = 0;
  return q;
}

void queue_destroy(Queue *q) {
  assert(q != NULL || q->nodes != NULL);
  free(q->nodes);
  free(q);
  q = NULL;
}

void queue_add(Queue *q, int e) {
  assert( !queue_full(q) );
  q->nodes[q->end] = e;
  q->end++;
  q->end %= q->size;
}

int queue_get(Queue *q) {
  int r = q->nodes[q->begin];
  assert( !queue_empty(q) );
  r = q->nodes[q->begin];
  q->begin++;
  q->begin %= q->size;
  return r;
}

// int main() { return 0; }
