typedef struct bNode{
    unsigned long index;
    int state;
    int neibor_num;
    int married;
    unsigned long  neighbors[100];
    unsigned long  partner;
}bNode;
typedef struct tNode{
    unsigned long index;
    int state;
    int neibor_num;
    int married;
    unsigned long  neighbors[100];
    unsigned long  partner;
}tNode;

bNode *bottoms,*xbottoms;
tNode *tops,*xtops;

typedef struct TimeNode{
    int events; // 0 recover 1 infect
    int bort; // b 0 t 1
    double t; // events happen time
    struct TimeNode *next; //next node
    unsigned long rm_node;
    unsigned long inf_node;
}TimeNode;
TimeNode *Fpointer;
void initFpointer(){
    Fpointer=malloc(sizeof(TimeNode));
    Fpointer->next=NULL;
}
void AddTimeRm(unsigned long inf_n, double rec_t,int bt){
    
    TimeNode *Node;
    Node=malloc(sizeof(TimeNode));
    Node->events=0;
    Node->bort=bt;
    Node->t=rec_t;
    Node->next=NULL;
    Node->rm_node=inf_n;
    Node->inf_node=-1;
    if (Fpointer->next==NULL) Fpointer->next=Node;
    else{
        TimeNode *tmp=Fpointer->next;
        if (tmp->next==NULL){
            if (Node->t>tmp->t) tmp->next=Node;
            else {
                Node->next=tmp;
                Fpointer->next=Node;
            }
        }
        else if(Node->t<=tmp->t){
            Node->next=tmp;
            Fpointer->next=Node;
        }
        else{
            int flag=0;
            while(tmp->next!=NULL){
                if(Node->t<tmp->next->t){
                    Node->next=tmp->next;
                    tmp->next=Node;
                    flag=1;
                    break;
                }
                tmp=tmp->next;
            }
            if (flag==0){
                tmp->next=Node;
            }
        }

    }
    
}
void AddTimeInf(unsigned long inf_n,double inf_t,int bt){

    TimeNode *Node;
    Node=malloc(sizeof(TimeNode));
    Node->events=1;
    Node->bort=bt;
    Node->t=inf_t;
    Node->next=NULL;
    Node->rm_node=-1;
    Node->inf_node=inf_n;
    if (Fpointer->next==NULL) Fpointer->next=Node;
    else{
        TimeNode *tmp=Fpointer->next;
        if (tmp->next==NULL){
            if (Node->t>tmp->t) tmp->next=Node;
            else {
                Node->next=tmp;
                Fpointer->next=Node;
            }
        }
        else if(Node->t<=tmp->t){
            Node->next=tmp;
            Fpointer->next=Node;
        }
        else{
            int flag=0;
            while(tmp->next!=NULL){
                if(Node->t<tmp->next->t){
                    Node->next=tmp->next;
                    tmp->next=Node;
                    flag=1;
                    break;
                }
                tmp=tmp->next;
            }
            if (flag==0){
                tmp->next=Node;
            }
        }

    }
}
TimeNode* TakeMin(){
    if (Fpointer->next!=NULL) return Fpointer->next;
    else {
        printf("The Min node Pointer is NULL!\n");
    }
}
void RemoveNodeMin(){
    TimeNode *tmp=Fpointer->next;
    Fpointer->next=tmp->next;
    free(tmp);
}
void Infbottoms(unsigned long inf_n,double t){
    double rec_t=gsl_ran_exponential(r,1.0/gamm);
    AddTimeRm(inf_n,rec_t+t,0);
    int j;
    for(j=0;j<bottoms[inf_n].neibor_num;j++){
        if (tops[bottoms[inf_n].neighbors[j]].state==0){
            if (bottoms[inf_n].married==1){
                if (bottoms[inf_n].neighbors[j]==bottoms[inf_n].partner){
                    if (gsl_ran_flat(r,0,1)<(1-exp(-beta01*rho*rec_t)))
                    AddTimeInf(bottoms[inf_n].partner,gsl_ran_exponential(r,1.0/(beta01*rho))+t,1);
                } 
                else {
                    if(gsl_ran_flat(r,0,1)<(1-exp(-beta01*alpha*rec_t)))
                    AddTimeInf(bottoms[inf_n].neighbors[j],gsl_ran_exponential(r,1.0/(beta01*alpha))+t,1);
                }
            }
            else{
                if(gsl_ran_flat(r,0,1)<(1-exp(-beta01*rec_t)))
                AddTimeInf(bottoms[inf_n].neighbors[j],gsl_ran_exponential(r,1.0/beta01)+t,1);
            }
        }
    }
}
void Inftops(unsigned long inf_n,double t){
    double rec_t=gsl_ran_exponential(r,1.0/gamm);
    AddTimeRm(inf_n,rec_t+t,1);
    int j;
    for(j=0;j<tops[inf_n].neibor_num;j++){
        if (bottoms[tops[inf_n].neighbors[j]].state==0){
            if(tops[inf_n].married==1){
                 if (tops[inf_n].neighbors[j]==tops[inf_n].partner){
                    if (gsl_ran_flat(r,0,1)<(1-exp(-beta10*rho*rec_t)))
                    AddTimeInf(tops[inf_n].partner,gsl_ran_exponential(r,1.0/(beta10*rho))+t,0);
                } 
                else {
                    if(gsl_ran_flat(r,0,1)<(1-exp(-beta10*alpha*rec_t)))
                    AddTimeInf(tops[inf_n].neighbors[j],gsl_ran_exponential(r,1.0/(beta10*alpha))+t,0);
                }
            }
            else{
                if(gsl_ran_flat(r,0,1)<(1-exp(-beta10*rec_t)))
                AddTimeInf(tops[inf_n].neighbors[j],gsl_ran_exponential(r,1.0/beta10)+t,0);
            }
        }
       
    }
}

