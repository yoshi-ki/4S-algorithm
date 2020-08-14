#include "stdio.h"
#include "stdlib.h"

struct Stack{
  int position;
  int* array;
  unsigned int capacity;
};

struct Stack* create(unsigned int capacity){
  struct Stack* stack = (struct Stack*)malloc(sizeof(struct Stack));
  stack -> position = -1;
  stack -> capacity = capacity;
  stack -> array = (int*)malloc(stack->capacity * sizeof(int));
  return stack;
};

void push(int v, struct Stack* s){
  s->array[++s->position] = v;
  return;
};

int pop(struct Stack* s){
  return s->array[--s->position];
};



int main(){
  struct Stack* s = create(100000);
  push(1,s);
  push(2,s);
  push(3,s);
  printf("%d\n",pop(s));
  printf("%d\n",pop(s));

  return 0;
}
