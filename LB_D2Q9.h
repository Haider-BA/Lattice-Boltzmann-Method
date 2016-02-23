#ifndef LB_D2Q9_H
#define LB_D2Q9_H
#endif

#include<math.h>
#include<glib.h>
#include<stdio.h>
#include"NACA00xx.h"

typedef struct NODE{
    //int num;
    int i; //grid
    int j;
    double q[9]; //directoin (x,y) -> (i,j); (x,y)=(i,j)+t*v[][n];  
    int s; // count
}Node;

void equilibrium(double feq[9],double rho, double ux,double uy){
    int n=0;
    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};
    double W[9] = {4./9,
                   1./9,  1./9,  1./9,  1./9,
                   1./36, 1./36, 1./36, 1./36,};

    for(n=0;n<9;n++){
        feq[n] = rho*W[n]*(2-sqrt(1+3*ux*ux))*(2-sqrt(1+3*uy*uy))*pow((2*ux+sqrt(1+3*ux*ux))/(1-ux),v[0][n])*pow((2*uy+sqrt(1+3*uy*uy))/(1-uy),v[1][n]); 
    }
    //printf("%f ",rho);
    //printf("\n");

}

void collide(int Nx, int Ny, double f[Nx][Ny][9], double alpha, double beta){ 
    double rho = 0;
    double jx = 0;
    double jy = 0;
    double ux = 0;
    double uy = 0;
    int n;
    int i,j;

    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};
    double W[9] = {4./9,
                   1./9,  1./9,  1./9,  1./9,
                   1./36, 1./36, 1./36, 1./36,};
    double feq;
    
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            
            for(n=0;n<9;n++){
                rho += f[i][j][n];
            }
                      
            for(n=0;n<9;n++){
                jx += f[i][j][n]*v[0][n];
                jy += f[i][j][n]*v[1][n];

            }
            ux = jx/rho;
            uy = jy/rho;
            for(n=0;n<9;n++){
                feq = rho*W[n]*(2-sqrt(1+3*ux*ux))*(2-sqrt(1+3*uy*uy))*pow((2*ux+sqrt(1+3*ux*ux))/(1-ux),v[0][n])*pow((2*uy+sqrt(1+3*uy*uy))/(1-uy),v[1][n]); 
                f[i][j][n] = f[i][j][n] + alpha*beta*(feq-f[i][j][n]);
            }
            //printf("%f %f %f %f \n",rho,ux,uy,sqrt(ux*ux+uy*uy)*sqrt(3.0));
            rho = 0;
            jx = 0;
            jy = 0;
        }
    }

}


void advect(int Nx, int Ny, double fdt[Nx][Ny][9],double f[Nx][Ny][9]){
    int n;
    int i,j;
    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};

    // copy 
    /*for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            for(n=0;n<9;n++){
                ft[i][j][n] = f[i][j][n];
            }
        }
    }*/
    int idt,jdt;
    // advect without b.c.
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
                    for(n=0;n<9;n++){
                        idt = i+v[0][n];
                        jdt = j+v[1][n];
                        if((idt>-1)&&(idt<Nx)&&(jdt>-1)&&(jdt<Ny)){
                            fdt[idt][jdt][n] = f[i][j][n];
                        }
                    }
        }
    }

}


void bc_periodic(int Nx, int Ny, double fdt[Nx][Ny][9], double f[Nx][Ny][9], GList* list){
    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};
    int n,i,j;
    int ri,rj;
    GList* iterator = NULL;
    
    for(iterator = list; iterator; iterator = iterator->next){
        i = ((Node *)iterator->data)->i;
        j = ((Node *)iterator->data)->j;
        //printf("bc points: %d %d ",i,j);
        for(n=0;n<9;n++){
            ri = (i+v[0][n]+Nx)%Nx;
            rj = (j+v[1][n]+Ny)%Ny;
            fdt[ri][rj][n] = f[i][j][n];
            //printf(" %d %d",ri,rj);
        }
        //printf("\n");
    }

} 

void bc_tube(int Nx, int Ny, double fdt[Nx][Ny][9], double f[Nx][Ny][9],double ux0, double uy0, double rho0){
    //int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 //{0, 0, 1,  0, -1, 1,  1, -1, -1},};
                //0, 1, 2,  3,  4, 5,  6,  7, 8
    //int vsx[Q] = {0, 1, 4, 3, 2, 8, 7, 6, 5}; // y->-y
    //int vsx_fw[Q] = {0, 1, 2, 3, 2, 5, 6, 6, 5}; // y->|y| i.e. +y->+y -y->+y
    //int vsx_wf[Q] = {0, 1, 4, 3, 4, 8, 7, 7, 8}; // y->-|y| i.e.+y->-y -y->-y
    double feq[9];
    int i,j;

    // bottom: y=0 free-slip: -y->+y;
    j=0;
    for(i=0;i<Nx;i++){
        fdt[i][j][2] = f[i][j][4];
        fdt[i][j][5] = f[i][j][8];
        fdt[i][j][6] = f[i][j][7];
    }
    // top: y=Ny-1 free-slip: +y->-y;
    j=Ny-1;
    for(i=0;i<Nx;i++){
        fdt[i][j][4] = f[i][j][2];
        fdt[i][j][8] = f[i][j][5];
        fdt[i][j][7] = f[i][j][6];
    }
    // left: x=0 equilibrium: +x
    //equilibrium(feq,rho0,ux0,0);
    i=0;
    for(j=0;j<Ny;j++){
        
        equilibrium(feq,rho0,ux0,uy0);
        fdt[i][j][1] = feq[1];
        fdt[i][j][5] = feq[5];
        fdt[i][j][8] = feq[8];
    }
    // right: x=Nx-1 last step
    i=Nx-1;
    for(j=0;j<Ny;j++){
        fdt[i][j][6] = f[i][j][6];
        fdt[i][j][3] = f[i][j][3];
        fdt[i][j][7] = f[i][j][7];
    }

}


void bc_grad(int Nx, int Ny, double fdt[Nx][Ny][9], double f[Nx][Ny][9], GList* list, Airfoil NACA00xx, double rho0, double beta){
    int i,j,k;
    int m,n;
    int s;
    double q;
    

    double ux,uy;
    double rho=0;
    double jx=0;
    double jy=0;
    double cs=1./sqrt(3.0);

    //double rho_s[9];  //surrounding density
    double u_s[9][2]; //surrounding velocity
    double dudx_s[2][2]; //surrounding gradient

    double rho_t; // target density
    double ux_t; // target x velocity
    double uy_t; // targe y velocity
    double dudx_t[2][2]; // target gradient

    double Pneq_t[2][2]; 
    double Peq_t[2][2];
    double P_t[2][2];

    double xc=NACA00xx.xc;
    double yc=NACA00xx.yc;
    double theta=NACA00xx.theta;
    double ux_c=NACA00xx.ux;
    double uy_c=NACA00xx.uy;
    double wz_c=NACA00xx.wz;
    double r;

    double ux_w; // wall point velocity
    double uy_w; 
    
    //double rho0=1.0; // initial density

    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};
                //0, 1, 2,  3,  4, 5,  6,  7, 8
                //0, x, y, -x, -y, xy, -xy,-x-y,x-y
    double W[9] = {4./9,
                   1./9,  1./9,  1./9,  1./9,
                   1./36, 1./36, 1./36, 1./36,};

    int vbb[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};


    GList* iterator;


    for(iterator=list;iterator;iterator=iterator->next){
        i = ((Node*)iterator->data)->i;
        j = ((Node*)iterator->data)->j;
        s = ((Node*)iterator->data)->s;
         
        for(m=0;m<9;m++){
            for(n=0;n<9;n++){
                rho += f[i+v[0][m]][j+v[1][m]][n];
                jx += f[i+v[0][m]][j+v[1][m]][n]*v[0][n];
                jy += f[i+v[0][m]][j+v[1][m]][n]*v[1][n];
            }
            ux = jx/rho;
            uy = jy/rho;
            
            //rho_s[m] = rho;
            u_s[m][0] = ux;
            u_s[m][1] = uy;

            rho=0;
            jx=0;
            jy=0;
        }

        dudx_s[0][0] = (u_s[1][0]-u_s[3][0])/2;  //dux / dx
        dudx_s[1][0] = (u_s[1][1]-u_s[3][1])/2;  //duy / dx
        dudx_s[0][1] = (u_s[2][0]-u_s[4][0])/2;  //dux / dy
        dudx_s[1][1] = (u_s[2][1]-u_s[4][1])/2;  //duy / dy

        dudx_t[0][0] = dudx_s[0][0];
        dudx_t[1][0] = dudx_s[1][0];
        dudx_t[0][1] = dudx_s[0][1];
        dudx_t[1][1] = dudx_s[1][1];
        
        rho_t = 0;
        ux_t = 0;
        uy_t = 0;

        for(n=0;n<9;n++){
            q = ((Node*)iterator->data)->q[n];
            if((q>0)&&(q<1)){
                r = sqrt(pow(i+q*v[0][n]-xc,2)+pow(j+q*v[1][n]-yc,2));
                ux_w = ux_c + wz_c*r*-sin(theta);
                uy_w = uy_c + wz_c*r*cos(theta);

                k=vbb[n]; // point from solid to fluid
                rho_t += f[i][j][n];// rho bounce back
                rho_t += 6*W[k]*rho0*(v[0][k]*ux_w+v[1][k]*uy_w); // rho local
                ux_t += 1./s*(q*u_s[k][0]+ux_w)/(1+q);
                uy_t += 1./s*(q*u_s[k][1]+uy_w)/(1+q);

                if((n==1)||(n==3)){ // x direction intersection
                    dudx_t[0][0] = (u_s[k][0]-ux_w)/(v[0][k]*(1+q));
                    dudx_t[1][0] = (u_s[k][1]-uy_w)/(v[0][k]*(1+q));
                }
                if((n==2)||(n==4)){
                    dudx_t[0][1] = (u_s[k][0]-ux_w)/(v[1][k]*(1+q));
                    dudx_t[1][1] = (u_s[k][1]-uy_w)/(v[1][k]*(1+q));
                }


            }else{
                k=vbb[n];
                rho_t += f[i][j][k];
            }
        }
        printf("%d %d %d %f %f %f\n",i,j,s,rho_t,ux_t,uy_t);

        Peq_t[0][0] = rho_t*pow(cs,2) + rho_t*ux_t*ux_t; // x/x
        Peq_t[1][0] = rho_t*uy_t*ux_t;
        Peq_t[0][1] = rho_t*ux_t*uy_t;
        Peq_t[1][1] = rho_t*pow(cs,2) + rho_t*uy_t*uy_t;

        Pneq_t[0][0] = -rho_t*pow(cs,2)/(2*beta)*(dudx_t[0][0] + dudx_t[0][0]);
        Pneq_t[1][0] = -rho_t*pow(cs,2)/(2*beta)*(dudx_t[1][0] + dudx_t[0][1]);
        Pneq_t[0][1] = -rho_t*pow(cs,2)/(2*beta)*(dudx_t[0][1] + dudx_t[1][0]);
        Pneq_t[1][1] = -rho_t*pow(cs,2)/(2*beta)*(dudx_t[1][1] + dudx_t[1][1]);
        P_t[0][0] = Peq_t[0][0] + Pneq_t[0][0];
        P_t[1][0] = Peq_t[1][0] + Pneq_t[1][0];
        P_t[0][1] = Peq_t[0][1] + Pneq_t[0][1];
        P_t[1][1] = Peq_t[1][1] + Pneq_t[1][1];

        for(n=0;n<9;n++){
            q = ((Node*)iterator->data)->q[n];
            if((q>0)&&(q<1)){
                k=vbb[n];
                fdt[i][j][k] = W[k]*(rho_t + rho_t/pow(cs,2)*(ux_t*v[0][k]+uy_t*v[1][k]) + 1/(2*pow(cs,4))*((P_t[0][0]-rho_t*pow(cs,2))*(v[0][k]*v[0][k]-pow(cs,2))+P_t[1][0]*v[1][k]*v[0][k]+P_t[0][1]*v[0][k]*v[1][k]+(P_t[1][1]-rho_t*pow(cs,2))*(v[1][k]*v[1][k]-pow(cs,2))));
            }
        }

    }
}



  






   





