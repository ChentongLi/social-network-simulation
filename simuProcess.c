#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define BNodeNum 10000
#define TNodeNum 8000
#define mu 10

gsl_rng *r;
//0 and 1 have different infection rate 9:1 (?variable)
//1 have a probability to connect a woman (10 millions)
//overlap 1 000 1
//1-0 soulmate decrease the connect rate
//social network development
//thailand connected rate decrease
//rate from 00 01 11
//circumsed rate
double marriedrate=0.6;
double beta10=0.04;
double beta01=0.02;
double gamm=0.03;
//double gam=2.2;
double rho=2.0;
double alpha=0.2;
#include "TimeNode.h"

int isNeighbor(unsigned long *neighbor,int num,unsigned long key){
    int i;
    for (i=0;i<num;i++)
        if(neighbor[i]==key)
            return 1;
    return 0;
}
int BuildGraph(){

    unsigned long  i;
    int j;
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    xbottoms=calloc(BNodeNum,sizeof(bNode));
    xtops=calloc(TNodeNum,sizeof(tNode));

    unsigned long *bedge;
    unsigned long *tedge;
    unsigned long bL=mu*BNodeNum,tL=mu*TNodeNum;
    unsigned long bk=0,tk=0;
    bedge=malloc(bL*sizeof(unsigned long));
    tedge=malloc(tL*sizeof(unsigned long));

    for (i=0;i<BNodeNum;i++){
        xbottoms[i].index=i;
        xbottoms[i].state=0;
        xbottoms[i].neibor_num=0;
        int tmpk=gsl_ran_poisson(r,mu); //neibor num
        //int tmpdeg=pow(gsl_rng_uniform(r),-1/gam)+1;
        //int tmpk=tmpdeg>99? 99:tmpdeg;
        xbottoms[i].married=0;
        xbottoms[i].partner=-1;
        int j;
        for(j=0;j<tmpk;j++){
            if (bk>bL-1) {
                bL+=100;
                bedge=realloc(bedge,bL*sizeof(unsigned long));
            }
            bedge[bk]=i;
            bk++;
        }
    }
    for (i=0;i<TNodeNum;i++){
        xtops[i].index=i;
        xtops[i].state=0;
        xtops[i].neibor_num=0;
        int tmpk=gsl_ran_poisson(r,mu); //neibor num
        //int tmpdeg=pow(gsl_rng_uniform(r),-1/gam)+1;
        //int tmpk=tmpdeg>99? 99:tmpdeg;
        xtops[i].married=0;
        xtops[i].partner=-1;
        int j;
        for(j=0;j<tmpk;j++){
            if (tk>tL-1) {
                tL+=100;
                tedge=realloc(tedge,tL*sizeof(unsigned long));
            }
            tedge[tk]=i;
            tk++;
        }
    }
    int fail=0;
    while (tk>0 && bk>0 && fail<100){
        unsigned long bii=gsl_rng_get(r)%bk;
        unsigned long tii=gsl_rng_get(r)%tk;
        unsigned long ti=tedge[tii];
        unsigned long bi=bedge[bii];
        if (isNeighbor(xtops[ti].neighbors,xtops[ti].neibor_num,bi)) fail++;
        else {
            xtops[ti].neighbors[xtops[ti].neibor_num]=bi;
            xtops[ti].neibor_num++;
            xbottoms[bi].neighbors[xbottoms[bi].neibor_num]=ti;
            xbottoms[bi].neibor_num++;
            unsigned long kk;
            tk--;
            for(kk=tii;kk<tk;kk++) tedge[kk]=tedge[kk+1];
            bk--;
            for(kk=bii;kk<bk;kk++) bedge[kk]=bedge[kk+1];
            fail=0;
        }
    }
    free(tedge);
    free(bedge);

    for(i=0;i<TNodeNum;i++){
        if(gsl_ran_flat(r,0,1)>marriedrate) xtops[i].married=0;
        else xtops[i].married=1;
        if(xtops[i].neibor_num>0 && xtops[i].married==1){
            int ii=gsl_rng_get(r)%xtops[i].neibor_num;
            unsigned long partner=xtops[i].neighbors[ii];
            if(xbottoms[partner].married==0){
                xtops[i].partner=partner;
                xbottoms[partner].married=1;
                xbottoms[partner].partner=i;
            }
            else xtops[i].married=0;
        }
    }
   // printf("%lu\n",k);
    printf("\nnetwork builded!\n");
    return 1;
}
int copygraph(){
    unsigned long i;
    int j;
    bottoms=calloc(BNodeNum,sizeof(bNode));
    tops=calloc(TNodeNum,sizeof(tNode));

    for(i=0;i<TNodeNum;i++){
        tops[i].index=xtops[i].index;
        tops[i].state=xtops[i].state;
        tops[i].married=xtops[i].married;
        tops[i].neibor_num=xtops[i].neibor_num;
        tops[i].partner=xtops[i].partner;
        for(j=0;j<xtops[i].neibor_num;j++){
            tops[i].neighbors[j]=xtops[i].neighbors[j];
        }
    }
    for(i=0;i<BNodeNum;i++){
        bottoms[i].index=xbottoms[i].index;
        bottoms[i].state=xbottoms[i].state;
        bottoms[i].married=xbottoms[i].married;
        bottoms[i].neibor_num=xbottoms[i].neibor_num;
        bottoms[i].partner=xbottoms[i].partner;
        for(j=0;j<xbottoms[i].neibor_num;j++){
            bottoms[i].neighbors[j]=xbottoms[i].neighbors[j];
        }
    }
    return 1;
}
int destroygraph(){
    free(tops);
    free(bottoms);
    return 1;
}

void simulate(){

    unsigned long inf_bnum=5;
    unsigned long inf_tnum=5;
    double t=0.0;
    int n=0;
    initFpointer();
    for (n=0;n<inf_bnum;n++){
        unsigned long inf_n=gsl_rng_get(r)%BNodeNum;
        bottoms[inf_n].state=1;
        Infbottoms(inf_n,0.0);
    }
    for (n=0;n<inf_tnum;n++){
        unsigned long inf_n=gsl_rng_get(r)%TNodeNum;
        tops[inf_n].state=1;
        Inftops(inf_n,0.0);
    }
    FILE *fp;
    char name[20];
    //sprintf(name,"alpha_%.3lf.csv",alpha); //print out name
    sprintf(name,"rho_%.3lf.csv",rho);
    //sprintf(name,"beta10_%.3lf.csv",beta10);
    //sprintf(name,"beta01_%.3lf.csv",beta01);
    //sprintf(name,"gamma_%.3lf.csv",gamm);
    fp=fopen(name,"w+");
    fprintf(fp,"Total,Set1,Set2,time\n");
    int flag=0;
    while(t<1000){
        TimeNode *tmpNode=TakeMin();
        if(tmpNode->events==0){
            if (tmpNode->bort==0){
                if (bottoms[tmpNode->rm_node].state==1){
                    inf_bnum--;
                    bottoms[tmpNode->rm_node].state=2;
                }
                else flag=1;

            }
            else{
                if(tops[tmpNode->rm_node].state==1){
                    inf_tnum--;
                    tops[tmpNode->rm_node].state=2;
                }
                else flag=1;
            }
        }
        else{
            if(tmpNode->bort==0){
                if (bottoms[tmpNode->inf_node].state==0){
                    inf_bnum++;
                    bottoms[tmpNode->inf_node].state=1;
                    Infbottoms(tmpNode->inf_node,tmpNode->t);
                }
                else flag=1;
            }
            else{
                if(tops[tmpNode->inf_node].state==0){
                    inf_tnum++;
                    tops[tmpNode->inf_node].state=1;
                    Inftops(tmpNode->inf_node,tmpNode->t);
                }
                else flag=1;
            }
        }
        if(flag==0){
            t=tmpNode->t;
            fprintf(fp,"%lu,%lu,%lu,%lf\n",inf_bnum+inf_tnum,inf_bnum,inf_tnum,tmpNode->t);
        }
        RemoveNodeMin();
        if ((inf_bnum+inf_tnum)==0 || Fpointer->next==NULL) break;
        flag=0;
    }
   fclose(fp);
}
int main(){
    srand(time(0));
    if(BuildGraph()){
        int i;
        for(i=0;i<5;i++){
            //alpha=0.2*i+0.01;  //here is the place you should change
//alpha is the parameter about when two person get married the rate of one of them have sex with others will decrease
//rho is they have sex with each other will increase
//beta10 is the infect rate from 1 to 0
//gamm is the recover rate
            rho=0.25*i+1;
            //beta10=0.01*i+0.04;
            //beta01=0.003*i+0.02;
            //gamm=0.01*i+0.01;
            copygraph();
            simulate();
            destroygraph();
            printf("finished %dth \n",i);
        }
    }
    return 0;
}

