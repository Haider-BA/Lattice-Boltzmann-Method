#include<stdio.h>
#include<glib.h>
#include<stdbool.h>
//#include"NACA00xx.h"
//#include"NACA00xx_A.h"
#include"LB_D2Q9_A.h"
#define M_PI 3.14159265358979323846

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
    int s=0; // count
    bool flag=0; // flag=1; append the node
    //int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                // {0, 0, 1,  0, -1, 1,  1, -1, -1},};
    //int ap[Nx][Ny]; // point couting
    //double qw[Nx][Ny][9]; 
    /*for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            ap[i][j]=0;
            qw[i][j][0]=0;
            for(n=1;n<9;n++){
                qw[i][j][n]=1;
            }
        }
    }*/
    //printf("theta=%f\n",theta);
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
                 
                Node* node = g_new(Node,1);
		node->i = i;
		node->j = j;
                for(n=0;n<9;n++){
			if(0==n){
				node->q[n]=0;
			}else{
				node->q[n]=1;
			}
                    d = NR(point,theta,c,n,xx);
                    if((d>0)&&(d<1)){
			    s++;
			    node->q[n]=d;
			    flag=1;
		    }
                        //ap[i][j]++;
                        //qw[i][j][n]=d;
                        // global:xw,yw,xc,yc
                        //xw = gx + d*v[0][n];
                        //yw = gy + d*v[1][n];
                        //r = sqrt(pow(xw-xc,2)+pow(yw-yc,2));
                        //vw[i][j][n][0] = ux + wz*r*-sin(theta); //ki=>-j
                        //vw[i][j][n][1] = uy + wz*r*cos(theta);  //kj=>i
                }
		if(flag==1){
			node->s = s;
                	list = g_list_append(list,node);
                }else{
			g_free(node);
		}
		flag=0;
		s=0;

            }
        }
    }
    
    /*for(i=0;i<Nx;i++){
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
    }*/

    
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
    
    //printf("list length: %d\n", g_list_length(list));
    return list;

}

int main(int argc, char* argv[]){
    int Nx=500;
    int Ny=500;
    Airfoil NACA00xx;
    //NACA00xx.xc=200;
    //NACA00xx.yc=Ny/2;
    //NACA00xx.theta=0.0;
    //NACA00xx.ux=0;
    //NACA00xx.uy=0;
    //NACA00xx.wz=0;
    NACA00xx.xx=0.10;
    NACA00xx.c=40;   // chord
    NACA00xx.xm=NACA00xxCent(NACA00xx.c,NACA00xx.xx);
    NACA00xx.ym=0;
    //double S=NACA00xxInt(NACA00xx.c,NACA00xx.c,NACA00xx.xx); //wing area
    //printf("wing area: %f\n",S);
    int timesteps;
    int time;
    double f[Nx][Ny][9];
    double fdt[Nx][Ny][9];
    double alpha=2.0;
    double Re=75;
    double chord=NACA00xx.c;
    double u0=0.01;
    double ux0;
    double uy0;
    double rho0=1; // rho0 also set to 1 in bc_grad function 
    double visc;
    double beta; //cs2 = 1/3;

    GList* list=NULL;
    GList* iterator=NULL;
    GList* temp=NULL;
    FILE *filerho = fopen("NACA_RHO_F_test_2T.dat","w");
    FILE *fileux = fopen("NACA_UX_F_test_2T.dat","w");
    FILE *fileuy = fopen("NACA_UY_F_test_2T.dat","w");
    
    //FILE *filedl = fopen("NACA_F_test3_LD.dat","w");
    double F[2];
    //double Lt_max,Dt_max;
    double cd,cl;
    //double Lt,Dt;

    int i,j;
    int n;
    int s;
    double rho,jx,jy,ux,uy;
    double feq[9];
    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};
    //u0 = pi*f*A0
    double A0=2.8*chord;
    double T=M_PI*A0/u0; // 0.25Hz
    double theta0=M_PI/2;
    double beta0=M_PI/4;
    double phi=M_PI/4;
    
    //u0 = M_PI/T*A0; 
    
    visc=u0*chord/Re;
    beta=1/(6*visc+1); 
    timesteps=2*T;//10*chord/u0;

	ux0 = 0;
	uy0 = 0;
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            equilibrium(feq,rho0,ux0,uy0);
            for(n=0;n<9;n++){
                f[i][j][n]=feq[n];
            }
        }
    }

    
    for(time=0;time<timesteps;time++){
	NACA00xx.xc=Nx/2+A0/2*cos(2*M_PI/T*time);
    	NACA00xx.yc=Ny/2;
    	NACA00xx.theta=theta0+beta0*sin(2*M_PI/T*time+phi);
    	//printf("xc: %f yc: %f th: %f \n",NACA00xx.xc,NACA00xx.yc,NACA00xx.theta);
	NACA00xx.ux=-A0*M_PI/T*sin(2*M_PI/T*time);
    	NACA00xx.uy=0;
    	NACA00xx.wz=beta0*2*M_PI/T*cos(2*M_PI/T*time+phi);
    	//printf("ux: %f uy: %f wz: %f \n",NACA00xx.ux,NACA00xx.uy,NACA00xx.wz);


    	
	list = Solid(Nx,Ny,NACA00xx,list);
    	//printf("list length: %d\n", g_list_length(list));
	if(time>0){// at the beginning no solid to fluid node
		s2f_F(Nx,Ny,fdt,time-1,time,rho0,list,NACA00xx.c,NACA00xx.xx,NACA00xx.xm,NACA00xx.ym,A0,beta0,1./T,phi,Nx/2,Ny/2,theta0);
	}

        advect(Nx,Ny,fdt,f,list); // periodic b.c.
        //printf("%d advect\n",time);
        /*s=0;
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			for(n=0;n<9;n++){
				if(-1==fdt[i][j][n]){
					s++;
				}
			}
		}
	}
	printf("f -1 number: %d \n", s);*/
        //bc_tube(Nx,Ny,fdt,f,ux0,uy0,rho0);
        //bc_freeslip(Nx,Ny,fdt,f);
	
        bc_grad(Nx,Ny,fdt,f,list,NACA00xx,rho0,beta);
	//CD = F[0]/(0.5*rho0*pow(u0,2)*chord);
	//CL = F[1]/(0.5*rho0*pow(u0,2)*chord);
	cd = 1.4-cos(2*NACA00xx.theta);
	cl = 1.2*sin(2*NACA00xx.theta);
	//printf("%f %f %f %f %f %f %f %f %f\n",time*1./timesteps,F[0],F[1],Dt_max,Lt_max,F[0]/Dt_max/S, F[1]/Lt_max/S, F[0]/(0.5*u0*u0*S*cd),F[1]/(0.5*u0*u0*S*cl));
	//fprintf(filedl,"%d %f %f %f %f %f %f %f %f\n",time,F[0],F[1],Dt_max,Lt_max, F[0]/Dt_max/S, F[1]/Lt_max/S, cd,cl);

        //printf("%d boundary\n",time);	

        //force_val(Nx,Ny,fdt,list,F);
	//printf("%f %f %f %f %f\n",time*1./timesteps,F[0],F[1],cd,cl);


        collide(Nx,Ny,fdt,alpha,beta,list);
        //printf("%d collide\n",time);

        force_val(Nx,Ny,fdt,list,NACA00xx,F);
	//printf("%f %f %f %f %f\n",time*1./timesteps,F[0],F[1],cd,cl);

	//fprintf(filedl, "%d %f %f %f %f\n",time,F[0],F[1],cd,cl);
        
        for(i=0;i<Nx;i++){
            for(j=0;j<Ny;j++){
                for(n=0;n<9;n++){
                    f[i][j][n]=fdt[i][j][n];
                }
            }
        }
	if(time>0){
	iterator = g_list_first(list);
    	while(iterator){
		s = ((Node*)iterator->data)->s;
		temp = iterator->next;
		if(-1==s){
			list=g_list_remove_link(list,iterator);
			g_free(iterator->data);
			g_list_free(iterator);
		}
		iterator = temp;
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
 
