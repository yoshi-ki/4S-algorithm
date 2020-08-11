#include "stdio.h"
#include "stdlib.h"

struct Stack{
  int val;
  struct Stack* next;
};

struct Stack* create(){
  return (struct Stack*)malloc(sizeof(struct Stack));
};

void push(struct Stack **s, int data) {
    struct Stack *new = create();
    new->val = data;
    new->next = *s;
    *s = new;
}

// struct Stack* push(int v, struct Stack* s){
//   struct Stack* s1 = create();
//   s1->val = v;
//   s1->next = s;
//   return s1;
// };

int pop(struct Stack** s){
  int re = (*s)->val;
  struct Stack* s1 = (*s)->next;
  free(*s);
  *s = s1;
  return re;
};



int main(){
  struct Stack* s = create();
  push(&s,1);
  push(&s,2);
  push(&s,3);
  printf("%d\n",pop(&s));
  printf("%d\n",pop(&s));

  return 0;
}
