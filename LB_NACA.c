#include<stdio.h>
#include<glib.h>
//#include"NACA00xx.h"
#include"LB_D2Q9.h"

GList* Wall(int Nx, int Ny, Airfoil NACA00xx, GList* list){
    // airfoil
    double xx=NACA00xx.xx;
    double c=NACA00xx.c;
    double xc=NACA00xx.xc;
    double yc=NACA00xx.yc;
    double theta=NACA00xx.theta;
    //double ux=NACA00xx.ux;
    //double uy=NACA00xx.uy;
    //double wz=NACA00xx.wz;
    double xm=NACA00xx.xm;
    double ym=NACA00xx.ym;
    double t=xx*c;
// search the points around airfoil
// length: chord;
// width: thicknss;
//  ------------------------------------------------
//  -------------------------------------------------
    int dx=1;
    int dy=1;
    int n;
    int i,j; //global grid point
    double x,y; //local grid point
    double point[2];   
    double d;
    double gx,gy; //global grid point gx=i,gy=j if dx=1,dy=1
    //double r; // distance between wall point and mass center
    //int s=0; // count
    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};
    int ap[Nx][Ny]; // point couting
    double qw[Nx][Ny][9]; 
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            ap[i][j]=0;
            qw[i][j][0]=0;
            for(n=1;n<9;n++){
                qw[i][j][n]=1;
            }
        }
    }
    printf("theta=%f\n",theta);
    for(i=0;i<Nx;i++){
        gx=i*dx;
        for(j=0;j<Ny;j++){
            // origin:(xc,yc) theta
            // grid point:(i*dx,j*dy)
            gy = j*dy;
            x = (gx-xc)*cos(theta)+(gy-yc)*sin(theta);
            y = (gx-xc)*-sin(theta)+(gy-yc)*cos(theta);
            // local:x,y,point
            // contraint: x in [0,c]; y in [-t/2,t/2]
            // grid constraint: [-1,c+1] [-t/2-1,t/2+1]
            if((x+xm>=-1)&&(x+xm<=c+1)&&(fabs(y+ym)<=t/2.+1)){
                point[0] = x+xm;
                point[1] = y+ym;
                 
                for(n=0;n<9;n++){
                    d = NR(point,theta,c,n,xx);
                    if((d>0)&&(d<1)){
                        ap[i][j]++;
                        qw[i][j][n]=d;
                        // global:xw,yw,xc,yc
                        //xw = gx + d*v[0][n];
                        //yw = gy + d*v[1][n];
                        //r = sqrt(pow(xw-xc,2)+pow(yw-yc,2));
                        //vw[i][j][n][0] = ux + wz*r*-sin(theta); //ki=>-j
                        //vw[i][j][n][1] = uy + wz*r*cos(theta);  //kj=>i
                    }
                }
            }
        }
    }
    
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            if(ap[i][j]!=0){
                Node* node = g_new(Node,1);
                node->i = i;
                node->j = j; 
                for(n=0;n<9;n++){
                    node->q[n] = qw[i][j][n];
                } 
                node->s = ap[i][j];
                list = g_list_append(list,node);
            }
        }
    }

    
    /*FILE *filea = fopen("NACA_a.dat","w"); 
    FILE *filex = fopen("NACA_x.dat","w");
    FILE *filey = fopen("NACA_y.dat","w");
    
    GList* iterator=NULL;

    for(iterator=list;iterator;iterator=iterator->next){
        i = ((Node*)iterator->data)->i;
        j = ((Node*)iterator->data)->j;
        for(n=0;n<9;n++){
            d = ((Node*)iterator->data)->q[n];
            if((d>0)&&(d<1)){
                x = i + v[0][n]*d;
                y = j + v[1][n]*d;
                fprintf(filex,"%f \n",x);
                fprintf(filey,"%f \n",y);

            }
        }
        printf("%d %d",i,j);
        fprintf(filea,"%d %d\n",i,j);
    }*/
    
    printf("list length: %d\n", g_list_length(list));
    return list;

}

int main(int argc, char* argv[]){
    int Nx=200;
    int Ny=200;
    Airfoil NACA00xx;
    NACA00xx.xc=20;
    NACA00xx.yc=Ny/2;
    NACA00xx.theta=-0.2;
    NACA00xx.ux=0;
    NACA00xx.uy=0;
    NACA00xx.wz=0;
    NACA00xx.xx=0.12;
    NACA00xx.c=20;   // chord
    NACA00xx.xm=NACA00xxCent(NACA00xx.c,NACA00xx.xx);
    NACA00xx.ym=0;

    int timesteps;
    int time;
    double f[Nx][Ny][9];
    double fdt[Nx][Ny][9];
    double alpha=2.0;
    double Re=100;
    double chord=20;//NACA00xx.c;

    double u0=0.05;
    double rho0=1; // rho0 also set to 1 in bc_grad function 
    double visc;
    double beta; //cs2 = 1/3;

    GList* list=NULL;

    FILE *filerho = fopen("NACA_RHO.dat","w");
    FILE *fileux = fopen("NACA_UX.dat","w");
    FILE *fileuy = fopen("NACA_UY.dat","w");

    int i,j;
    int n;
    double rho,jx,jy,ux,uy;
    double feq[9];
    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};

    timesteps=10*chord/u0;
    visc=u0*chord/Re;
    beta=1/(6*visc+1); 
    list = Wall(Nx,Ny,NACA00xx,list);
    
    printf("list length: %d\n", g_list_length(list));

    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            equilibrium(feq,rho0,u0,0);
            for(n=0;n<9;n++){
                f[i][j][n]=feq[n];
            }
        }
    }
    
    for(time=0;time<timesteps;time++){
        advect(Nx,Ny,fdt,f);
        printf("%d advect\n",time);
        
        bc_tube(Nx,Ny,fdt,f,u0,rho0);
        bc_grad(Nx,Ny,fdt,f,list,NACA00xx,rho0,beta);
        printf("%d boundary\n",time);
        
        collide(Nx,Ny,fdt,alpha,beta);
        printf("%d collide\n",time);
        
        for(i=0;i<Nx;i++){
            for(j=0;j<Ny;j++){
                for(n=0;n<9;n++){
                    f[i][j][n]=fdt[i][j][n];
                }
            }
        }
    }
 
    rho = 0;
    jx = 0;
    jy = 0;
    ux = 0;
    uy = 0;
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
            //printf("%f %f %f",rho,ux,uy);
            fprintf(filerho,"%f ", rho);
            fprintf(fileux, "%f ", ux);
            fprintf(fileuy, "%f ", uy);
            rho = 0;
            jx = 0;
            jy = 0;

        
        }
        //printf("\n");
        fprintf(filerho, "\n");
        fprintf(fileux, "\n");
        fprintf(fileuy, "\n");
    }


    g_list_free(list);

    return 0;
}


/*int main(int argc, char* argv[]){
    int Nx=64;
    int Ny=64;
    Airfoil NACA0012;
    NACA0012.xc=20;
    NACA0012.yc=Ny/2;
    NACA0012.theta=0;
    NACA0012.ux=0;
    NACA0012.uy=0;
    NACA0012.wz=0;
    NACA0012.xx=0.12;
    NACA0012.c=20;   // chord
    NACA0012.xm=NACA00xxCent(NACA0012.c,NACA0012.xx);
    NACA0012.ym=0;

    GList* list=NULL;

    Wall(Nx,Ny,NACA0012,list);

    g_list_free(list);

    return 0;
}*/
 
