#ifndef LB_D2Q9_A_H
#define LB_D2Q9_A_H
#endif

#include<math.h>
#include<stdbool.h>
#include<glib.h>
#include<assert.h>
#include<stdio.h>
#include"NACA00xx_A.h"

// boundary node
typedef struct NODE{
    int i;
    int j;
    double q[9]; // distance
    int s; // number of intersection points. s=0: solid;s>0: boudary; s=-1: to be deleted;
}Node;

// search solid area including boundary points
GList* Solid(int Nx, int Ny, Airfoil NACA00xx, GList* list){
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
    int n;
    int i,j; //global grid point
    double x,y; //local grid point
    double point[2];   
    double d;
    //double r; // distance between wall point and mass center
    int s=0; // count
    bool b_flag=0; // flag=1; append the node (boundary)
    bool s_flag=0; // flag=1; append the node (solid)
    bool n_flag=0; // flag=1; solid to fluid nodes

    GList* iterator=NULL;
    GList* nlist=NULL; // new list: add new nodes
    GList* niterator=NULL; // new iterator
    GList* temp=NULL;
    int i_o,j_o; // old i,j;
    int s_o; // old s for solid to fluid node i.e. s_o=0 
    //printf("theta=%f\n",theta);
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // origin:(xc,yc) theta
            // grid point:(i*dx,j*dy)
            x = (i-xc)*cos(theta)+(j-yc)*sin(theta);
            y = (i-xc)*-sin(theta)+(j-yc)*cos(theta);
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
		}
                if(OBJ(0,point,theta,c,0,xx)<=0){
		for(n=0;n<9;n++){		   
                    d = NR(point,theta,c,n,xx); // NR function consider points outside the airfoil but inside the box.
                    if((d>0)&&(d<1)){
                        s++;
                        node->q[n]=d;
                        b_flag=1;
		    }
                }
		}
                if(OBJ(0,point,theta,c,0,xx)>0){
                    // OBJ>0: point locates inside the airfoil
                    // which is solid
                    s = 0;
                    s_flag=1;
                }

                if((1==b_flag)||(1==s_flag)){ // if point belongs to solid area includeing the boundary
                    node->s = s; // s>0 for boundary; s=0 for solid 
 		    assert(s>=0); 		    
                    nlist = g_list_prepend(nlist,node); // add new point

	    	    //printf("new: %d %d %d\n",i,j,s);

                }else{
                    g_free(node);
                }

                b_flag=0;
                s_flag=0;
                s=0;

            }
        }
    }
    iterator = g_list_first(list); 
    while(iterator){
	temp = iterator->next;
	i_o = ((Node*)iterator->data)->i;
        j_o = ((Node*)iterator->data)->j;
	s_o = ((Node*)iterator->data)->s;
	//printf("old: %d %d %d\n",i_o,j_o,s_o);
	assert(s_o>=0);
	    // check whether it is new point or old point
	    // nlist as a reference not change! 
	    // change list to leave only solid to fluid nodes
	    n_flag = 0;
	if(0==s_o){ // solid to fluid nodes
	    for(niterator=nlist;niterator;niterator=niterator->next){
		i = ((Node*)niterator->data)->i;
		j = ((Node*)niterator->data)->j;
		s = ((Node*)niterator->data)->s;	
		if((i==i_o)&&(j==j_o)){
			    n_flag = 1;
			    //printf("old: %d %d %d\n",i_o,j_o,s_o);
			    break;	
		}
		
	    }
	}
	if((1==n_flag)||(s_o>0)){
	    list = g_list_remove_link(list,iterator); // no longer belong to boundary points, deleted her
	    g_free(iterator->data);
	    g_list_free(iterator);
	    //printf("free data: boundary nodes\n");
	}else{  // keep the solid to fluid nodes	   
	    ((Node*)iterator->data)->s = -1; 
	}
	iterator = temp;
    }
    
    list = g_list_concat(list,nlist); // append new list   
    printf("list length: %d\n", g_list_length(list));
    return list;

}



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

void collide(int Nx, int Ny, double f[Nx][Ny][9], double alpha, double beta, GList* list){ 
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

    //GList* iterator=NULL;
    //double s;
    
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

    // collision doesn't happen in solid, set f to -1
    /*for(iterator=list;iterator;iterator=iterator->next){
        i = ((Node*)iterator->data)->i;
        j = ((Node*)iterator->data)->j;
        s = ((Node*)iterator->data)->s;
        if(0==s){ // 0:solid;-1:to be deleted nodes not change it here
            for(n=0;n<9;n++){
                f[i][j][n]=-1;
            }
        }
    }*/

}

void advect(int Nx, int Ny, double fdt[Nx][Ny][9],double f[Nx][Ny][9],GList* list){
    int n;
    int i,j;
    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};

    int idt,jdt;

    //GList* iterator=NULL;
    //int s;
    //double q;
    //int vbb[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    //int k;
    // advect fluid
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
                    for(n=0;n<9;n++){
                        idt = i+v[0][n];
                        jdt = j+v[1][n];
                        if((idt>-1)&&(idt<Nx)&&(jdt>-1)&&(jdt<Ny)){
                            fdt[idt][jdt][n] = f[i][j][n];
                        }
                        else{
                            // periodic b.c.
                            idt = (idt + Nx)%Nx;
                            jdt = (jdt + Ny)%Ny;
                            fdt[idt][jdt][n] = f[i][j][n];
                        }
                    }
            
        }
    }

    // not advect points in solid area
    // not advect some populatoins of boundary points (to be approximated by Grad)
    // set them to -1
    /*for(iterator=list;iterator;iterator=iterator->next){
        i = ((Node*)iterator->data)->i;
        j = ((Node*)iterator->data)->j;
        s = ((Node*)iterator->data)->s;
        if(0==s){ // 0: solid
            fdt[i][j][n] = -1;
        }
        if(s>0){ // >0: boundary
            for(n=0;n<9;n++){
                k = vbb[n];      
                q = ((Node*)iterator->data)->q[n];
                if((q>0)&&(q<1)){
                    fdt[i][j][k] = -1;
                }
            }
        }

    }*/

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


    GList* iterator=NULL;

    //Grad's approximation for boundary points and points from solid to fluid
    for(iterator=list;iterator;iterator=iterator->next){
        i = ((Node*)iterator->data)->i;
        j = ((Node*)iterator->data)->j;
        s = ((Node*)iterator->data)->s;
        if(s>0){   
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
        //printf("%d %d %d %f %f %f\n",i,j,s,rho_t,ux_t,uy_t);

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
        //if(-1==s){
        // old node from solid to fluid  
        //}

    }
}



void s2f_F(int Nx, int Ny, double f[Nx][Ny][9], double t0, double t1, double rho0, GList* list, double c, double xx, double xm, double ym, double A0, double B0, double freq, double phi, double x0, double y0, double alpha0){
// solid to fluid node for flapping wing
// assign equilibrium distributoin
    int i,j;
    int s;
    int n;

    double t=(t0+t1)*0.5;
    double uw[2];
    int grid[2];
    
    double feq[9];

    GList* iterator=NULL;
    for(iterator=list;iterator;iterator=iterator->next){
        i = ((Node*)iterator->data)->i;
        j = ((Node*)iterator->data)->j;
        s = ((Node*)iterator->data)->s;

        if(-1==s){
            assert((i>-1)&&(i<Nx));
            assert((j>-1)&&(j<Ny));

            grid[0] = i;
            grid[1] = j;
            Wall_F(t,t0,t1,grid,uw,c,xx,xm,ym,A0,B0,freq,phi,x0,y0,alpha0);

            equilibrium(feq,rho0,uw[0],uw[1]);

            for(n=0;n<9;n++){
              f[i][j][n] = feq[n];
            }
        }
    }

}

void force_val(int Nx, int Ny, double f[Nx][Ny][9], GList* list, Airfoil NACA00xx, double F[2]){
    int i,j,k,n;
    int s;
    double q;

    //force measure
    double jx_t;
    double jy_t;
    double cl; //cl=1.2*sin(2*theta)
    double cd; //cd=1.4-cos(2*theta)
    double Lt; //lift
    double Dt; //drag
    double Lmax=0;
    double Dmax=0;   
    
    double xc=NACA00xx.xc;
    double yc=NACA00xx.yc;
    double theta=NACA00xx.theta;
    double ux_c=NACA00xx.ux;
    double uy_c=NACA00xx.uy;
    double wz_c=NACA00xx.wz;
    double r=0;
    //printf("%f %f %f %f %f %f\n",xc,yc,theta,ux_c,uy_c,wz_c);
    double ux_w=0; // wall point velocity
    double uy_w=0; 
    int v[2][9]={{0, 1, 0, -1,  0, 1, -1, -1, 1},
                 {0, 0, 1,  0, -1, 1,  1, -1, -1},};
        
    int vbb[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    GList* iterator=NULL;
    //printf("list length: %d\n",g_list_length(list));

    cl = 1.2*sin(2*theta);
    cd = 1.4 - cos(2*theta);

    F[0] = 0;
    F[1] = 0;
    for(iterator=list;iterator;iterator=iterator->next){
        i = ((Node*)iterator->data)->i;
        j = ((Node*)iterator->data)->j;
        s = ((Node*)iterator->data)->s;
	//printf("%d %d %d \n",i,j,s);        
        //printf("%f %f %f\n",r,ux_w,uy_w);
        if(s>0){
        if((i<=1)||(i>=Nx-2)||(j<=1)||(j>=Ny-2)){
		printf("Airfoil touch the boundary, segmentation fault !");
	}
	
	jx_t = 0;
	jy_t = 0;

	for(n=0;n<9;n++){
            q = ((Node*)iterator->data)->q[n];
	    k = vbb[n];
	    //printf("q: %f\n",q); 
            if((q>0)&&(q<1)){
		r = sqrt(pow(i+q*v[0][n]-xc,2)+pow(j+q*v[1][n]-yc,2));
            	ux_w = ux_c + wz_c*r*-sin(theta);
            	uy_w = uy_c + wz_c*r*cos(theta);
		Dt = 0.5*(ux_w*ux_w+uy_w*uy_w)*cd;
		Lt = 0.5*(ux_w*ux_w+uy_w*uy_w)*cl;
		if(Dt>Dmax){
			Dmax = Dt;
		}
		if(Lt>Lmax){
			Lmax = Lt;
		}
		//printf("%d %d %d %f %f %f %f %f\n",i,j,s,r,ux_w,uy_w,Dt,Lt);
		jx_t += v[0][n]*(f[i][j][k]+f[i+v[0][n]][j+v[1][n]][n]);
    		jy_t += v[1][n]*(f[i][j][k]+f[i+v[0][n]][j+v[1][n]][n]);
	    }
	}
	//printf("jx: %f jy: %f\n",jx_t,jy_t);
	F[0] += jx_t;
	F[1] += jy_t;
	}

    }
	//F[0] = F[0]/Dmax;
	//F[1] = F[1]/Lmax;
}







